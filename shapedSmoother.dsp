declare name "shapedSmoother";
declare version "0.3.6";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

// AUC level-compensation multiplier (provides aucLevelMult, the dependency of
// aucLevelMultSwitched below). It is a degree-8 polynomial baked from
// aucLevelMultFormula — the source of truth, defined far below. Regenerate
// after any change to the curve or shapeMap with:
//   regen_auc_poly.sh shapedSmoother.dsp 8 257
import("auc_poly.lib");

// Commented out: only the reference LUT path (table_aucLevelMult, far below)
// needs aucTbl from this file. The live AUC multiplier uses auc_poly.lib above.
// import("auc_waveform.lib");

// ============================================================================
//  GUI — every user-facing control in one place.
//
//  Single source of truth for the group names: the three *Group functions
//  below. `vgroup("[1]Smoother", x)` is precisely what the old
//  "v:[1]Smoother/..." label prefix desugars to, and Faust merges groups by
//  path, so wrapping each control in SmootherGroup(...) yields an identical UI
//  tree. Renaming a group is now a one-line edit here.
//
//  Controls carry RAW slider values only; ms->s scaling and AUC compensation
//  live at the use sites (att / rel / rel_hold in the Parameters block). The
//  shape sliders (attackShapeSlider / releaseShapeSlider) keep their names so
//  the AUC level compensation (aucLevelMultFormula, baked into aucLevelMult in
//  auc_poly.lib) and env()'s shapeMap() can reference them unchanged.
// ============================================================================
MainGroup(x) = hgroup("[0]shapedSmoother", x);
TestGroup(x) = vgroup("[0]Test signal", x);
SmootherGroup(x) = vgroup("[1]Smoother", x);

// --- Test signal ---
testNoiseLevel = TestGroup(hslider("[0]noise level", 0, 0, 1, 0.001));
testNoiseRate = TestGroup(hslider("[1]noise rate", 42, 1, 1000, 1));
testBlockscale = TestGroup(hslider("[2]blockscale", 1, 0.01, 10, 0.01));
testFreq = TestGroup(hslider("[3]freq", 1, 0.001, 4, 0.001));
testStep1 = TestGroup(hslider("[4]step1", 0.75, -1, 1, 0.001));
testStep2 = TestGroup(hslider("[5]step2", 0.125, -1, 1, 0.001));

// --- Smoother (raw; scaled att/rel/rel_hold are derived in Parameters) ---
attMs = SmootherGroup(hslider("[0]att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01));
attackShapeSlider = SmootherGroup(hslider("[1]attack shape", 0, 0, 1, 0.001));
attAucComp = SmootherGroup(checkbox("[2]att auc comp"));
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));
relAucComp = SmootherGroup(checkbox("[5]rel auc comp"));
relHoldMs = SmootherGroup(hslider("[6]rel hold[unit:ms]", 50, 0, maxRelHold*1000, 1));

process = MainGroup((test2@(att_samples+brake_samples+rel_hold_samples), shapedSmoother(test2)));
shapedSmoother(x) = lookaheadX:env~(_, _, _)
    with {
        // Release hold, ported from lamb-rs (lamb.dsp `releaseHold` inside
        // `lookahead_compression_gain_mono`). Sits BEFORE the attack-lookahead
        // slidingMin, mirroring lamb's pipeline where releaseHold runs before
        // the SIN smoother's internal slidingMin.
        //
        // - `min(prevGain, rawX @ rel_hold_samples)` keeps the output
        //   monotonically non-rising: once we've gone down, we can't release.
        // - `slidingMin(rel_hold_samples+1, ...)` is the future-min within the
        //   hold window: we may go as low as that, but no lower is required.
        // - max() of those two yields the held value.
        //
        // Adds rel_hold_samples of latency on top of att_samples+brake_samples.
        held = x:(releaseHold~_);
        releaseHold(prevGain, rawX) = max(min(prevGain, rawX@rel_hold_samples),
            rawX:ba.slidingMin(rel_hold_samples+1, 1+maxRelHoldSamples));

        lookaheadX = held:ba.slidingMin(att_samples+1, 1+maxHoldSamples)@brake_samples;

        // brake: covers an additional brake_samples ahead, used purely for
        // early attack detection. Uses `held` so the brake sees the same
        // release-held signal as lookaheadX.
        lookahead_brake = held:ba.slidingMin(att_samples+brake_samples+1, 1+(maxHold+maxBrake)*maxSR);
        // ---- Hoisted per-shape constants (slider-rate) ----
        // Pure functions of the shape sliders. Faust's auto-lift doesn't push
        // these through `select2(releasing, sliderA, sliderB)`, so without
        // hoisting every `cheapCurveBase`/`curveScale`/`maxDerivativeBaseAttack`
        // call on `shape` would fire at audio rate. Precompute them per-shape.
        attackShape = shapeMap(attackShapeSlider);
        releaseShape = shapeMap(releaseShapeSlider);

        attackInvCurveScale = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackZero = cheapCurveBase(attackShape, 0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
        releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

        // Phase of the release derivative's peak (slider-rate). Used to pick the
        // accelerating (bottom) vs decelerating (top) branch robustly.
        releasePeakPhase = peakPhaseRelease(releaseShape);

        // v0.3.6: hoisted internals of cheapCurveBase for the cheapCurveBaseH
        // call in env() (slider-rate: pure functions of the shape sliders).
        attackCurveSqrtInner = sqrt((1/attackShape)-1);
        releaseCurveSqrtInner = sqrt((1/releaseShape)-1);
        attackCurveInvSqrtInner = 1/attackCurveSqrtInner;
        releaseCurveInvSqrtInner = 1/releaseCurveSqrtInner;
        attackCurveHalfInvC = 1/(2*attackShape);
        releaseCurveHalfInvC = 1/(2*releaseShape);

        // v0.3.6: per-sample 1/totalNRSteps hoisted to slider-rate. Bit-identical
        // (same operands and max(1) clamp, just computed once per block).
        attStep = 1/((att*ma.SR):max(1));
        relStep = 1/((rel*ma.SR):max(1));

        env(prev, prevPhase, prevTotalStep, lookaheadX) = result, newPhase, totalStep
            with {
                // Audio-rate scalars: select the precomputed value, no math.
                shape = select2(releasing, attackShape, releaseShape);
                // NB: arithmetic blend, NOT select2. `releasing` is 0/1, so this is
                // bit-identical to select2(releasing, attack, release) — but select2
                // traps its branches: Faust will not hoist a control-rate sub-term
                // across the select boundary, so `1/curveScale(shape)` kept firing an
                // atan() PER SAMPLE (the atan's arg was hoisted, the atan call wasn't).
                // As a `+ sel*(b-a)` blend both arms are ordinary operands, so the
                // whole reciprocal (atan + divide) lifts to slider-rate. Drops 2 atan
                // and 2 divides/sample. The sibling picks below don't carry a trapped
                // transcendental, so they stay as readable select2.
                invCurveScale = attackInvCurveScale+releasing*(releaseInvCurveScale-attackInvCurveScale);
                zeroVal = select2(releasing, attackZero, releaseZero);
                maxDerivVal = select2(releasing, attackMaxDerivBase, releaseMaxDerivBase);
                curveSqrtInner = select2(releasing, attackCurveSqrtInner, releaseCurveSqrtInner);
                curveInvSqrtInner = select2(releasing, attackCurveInvSqrtInner, releaseCurveInvSqrtInner);
                curveHalfInvC = select2(releasing, attackCurveHalfInvC, releaseCurveHalfInvC);

                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                step = select2(releasing, attStep, relStep);

                // make delta dependent on how much GR we have to do until the full lookahead, so it self corrects

                todo = lookaheadX-prev;

                // curve fraction covered at prevPhase (same form as gonnaDo's cb)
                cbPrev = (cheapCurveBaseH(shape, curveSqrtInner, curveInvSqrtInner, curveHalfInvC, select2(releasing, 1-prevPhase, prevPhase))-zeroVal)*invCurveScale;
                fracDone = select2(releasing, 1-cbPrev, cbPrev);

                // True the sample the aim point moves. lookaheadX is the value the env
                // chases, and slidingMin holds it bit-exactly constant while the window
                // minimum is unchanged, so `!= mem` fires precisely on a genuine target
                // change and is silent during a hold. (mem of an input is a 1-sample
                // delay; it does not enter the env feedback.)
                targetChanged = lookaheadX!=lookaheadX';

                // Span = (target - releaseStart). prev sits at fraction fracDone of the
                // running span, so releaseStart = prev - prevTotalStep*fracDone and the
                // live span is (target - prev) + prevTotalStep*fracDone = rawSpan.
                rawSpan = (lookaheadX-prev)+prevTotalStep*fracDone;
                // Attack keeps the old min-clamp; release uses the new target-change gate.
                //   attack : min(rawSpan, prevTotalStep) pins the span to prevTotalStep
                //     while rawSpan >= prevTotalStep (prevTotalStep/totalStep == 1, the
                //     resync round-trips, phase advances by +step => shaped curve
                //     preserved); bit-identical to the old line for a static target.
                //   release: pin to prevTotalStep while the target holds (same bit-exact
                //     resync), and adopt rawSpan only when the target actually moves. The
                //     old release branch used max(rawSpan, prevTotalStep), which pinned
                //     the static case the same way but could ONLY let the span grow, so
                //     once the release target moved closer the env kept aiming at the old,
                //     farther target and the transition couldn't shorten (the GR-goes-
                //     stationary bug). Keying off targetChanged lets the release span move
                //     either way.

                totalStep = select2(releasing,
                    ((lookaheadX-prev)+prevTotalStep*fracDone):min(prevTotalStep),
                    select2(targetChanged, prevTotalStep, rawSpan))*active;

                prevSpeed = prevTotalStep*derivativeBaseRelease(shape,
                    select2(releasing, 1-prevPhase, prevPhase));

                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0,
                    min(speedRatio,
                        maxDerivVal));

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                Readable_gonnaDo(phase) = (1-select2(releasing,
                    cheapCurveAttackS(shape, invCurveScale, zeroVal, phase),
                    cheapCurveReleaseS(shape, invCurveScale, zeroVal, phase)))*totalStep;
                // select2 evaluates BOTH arms; here they call cheapCurveBase at different
                // args (1-phase vs phase) so CSE can't merge them -> 2 transcendental
                // clusters/sample. Attack curve == 1 - release curve(1-x), so pick the
                // argument first and evaluate cheapCurveBase ONCE.
                gonnaDo(phase) = select2(releasing, cb, 1-cb)*totalStep
                    with {
                        cb = (cheapCurveBaseH(shape, curveSqrtInner, curveInvSqrtInner, curveHalfInvC, select2(releasing, 1-phase, phase))-zeroVal)*invCurveScale;
                    };

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio)))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                // Branch selection for `phaseAtMatchingSpeed`. ATTACK uses
                // `gonnaMakeIt` directly. RELEASE instead picks the accelerating
                // (bottom) inverse when either:
                //   prevPhase < releasePeakPhase  -> still before the derivative peak
                //   totalStep  > prevTotalStep    -> a larger target appeared (catch-up)
                // and the decelerating (top) inverse otherwise. Direct tests are
                // used here because for a single release step `projected` equals
                // the target, so `gonnaMakeIt` would be decided by rounding;
                // testing the conditions directly gives its intended value
                // without the tie. The phase resync (`phaseAtMatchingSpeed +
                // step`) is identical for both branches.
                releaseUseBottom = (prevPhase<releasePeakPhase)|(totalStep>prevTotalStep);
                phaseAtMatchingSpeed = select2(releasing,
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    select2(releaseUseBottom,
                        inverseDerivativeTopRelease(shape, clampedRatio),
                        inverseDerivativeBottomRelease(shape, clampedRatio)));

                speed = totalStep*derivativeReleaseS(shape,
                    invCurveScale,
                    select2(releasing, 1-newPhase, newPhase));

                braking = (prev>lookahead_brake)&(prev<=lookaheadX);
                brakeRamp = ((_+(1/brake_samples))*braking)~_:min(1);
                brakeMult = 1-brakeRamp*checkbox("enable brake");
                delta = speed*step*brakeMult;

                result = min(prev+delta, x@(att_samples+brake_samples+rel_hold_samples));
            };
    };

// Parameters
maxHold = 0.05;
maxBrake = 0.05;
maxRelHold = 0.05;
maxSR = 48000;

maxHoldSamples = maxHold*maxSR;
maxBrakeSamples = maxBrake*maxSR;
maxRelHoldSamples = maxRelHold*maxSR;
maxTotalSamples = (maxHold+maxBrake)*maxSR;

// AUC level compensation.
// Goal: changing the shape must NOT change the average level. A single-step
// response of duration T has time-area  T * AUC_attack(shape) * step. We scale
// the duration by  aucLevelMult(shape) = AUC_attack(shapeMap(1)) / AUC_attack(shape),
// so  T_eff * AUC = const  -> the area is shape-independent. The factor is
// normalized to <=1 (==1 at the sharpest shape), so durations only ever shorten.
// (aucLevelMult is the polynomial from auc_poly.lib; the exact factor it
// approximates is aucLevelMultFormula, defined far below.)
//
// We fold the factor straight into att/rel, so the lookahead, brake sizing AND
// fade speed all use the same compensated duration -> the single-step response
// stays a pure scaled curve (no early plateau). CONSEQUENCE: the lookahead, and
// therefore the reported plugin latency, now varies with shape (it is largest
// at the sharpest shape). The release-hold latency (rel_hold_samples) is a
// separate, shape-independent term added on top. If you need a FIXED latency
// instead, leave att/rel here un-scaled and instead multiply only `activeTime`
// inside env() by
// select2(releasing, aucLevelMult(attackShapeSlider), aucLevelMult(releaseShapeSlider));
// that keeps latency constant but reintroduces a shape-dependent plateau that
// only partially compensates the level.

// On/off for the AUC level compensation, independently for att and rel.
// The checkboxes (attAucComp / relAucComp) are declared in the GUI block.
// When off the multiplier is exactly 1.0, i.e. att/rel revert to the raw
// slider values — and the reported latency goes back to being
// shape-independent (see the AUC comment on the att/rel sliders).
//
// Blend the multiplier toward 1.0 when the switch is off. Branchless on
// purpose: `on` is 0 or 1 so it's exact, the multiplier is clamped to [0,1]
// so `on*(m-1)` can't blow up, and this is slider-rate regardless.
// aucLevelMult comes from auc_poly.lib (imported at the top).
aucLevelMultSwitched(on, s) = 1+on*(aucLevelMult(s)-1);

att = attMs/1000*aucLevelMultSwitched(attAucComp, attackShapeSlider);
att_samples = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
rel = relMs/1000*aucLevelMultSwitched(relAucComp, releaseShapeSlider);
rel_samples = rel*ma.SR:max(1);

// Release hold (ported from lamb-rs): linear scale, allows 0 (off).
// Default 50 ms matches lamb's `[08]release hold`.
rel_hold = relHoldMs*0.001;
rel_hold_samples = rel_hold*ma.SR;

// Brake window length per spec:
//   = 50 ms by default,
//   = max(att, rel) if both att AND rel are shorter than 50 ms.
// (min(50ms, max(att,rel)) gives exactly this.)
brake_samples = min(maxBrakeSamples, max(att_samples, rel_samples)):max(1);

test3 = os.lf_squarewave(testFreq);

// Test signal
test2 = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };
testSig = os.lf_sawpos(0.5)<:(((_>0.25)*testStep1)+((_>0.5)*testStep2));

// *************************************** the NEW curves: ******************************
// New curves by nuchi:
// https://www.wolframalpha.com/input?i=antiderivative+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29
// https://www.wolframalpha.com/input?i=inverse+function+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29+%2F+C
// https://www.desmos.com/calculator/9wtfhymvr0
// with scaling for c:
// https://www.desmos.com/calculator/dynyjjkuli
// with flip horizontal:
// https://www.desmos.com/calculator/bsr8cdn21v
cheapCurveBase(c, x) = (log(c*(x*x-1)+1)+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);
curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);
cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x) = cheapCurveRelease(c, x*-1+1)*-1+1;

// Variants taking precomputed 1/curveScale and cheapCurveBase(c,0), so they
// don't recompute curveScale and cheapCurveBase(c,0) per sample. Same result
// as cheapCurveRelease/Attack.
//   invScale = 1/curveScale(c)
//   zero     = cheapCurveBase(c, 0)
// v0.3.6: hoisted-term variant of cheapCurveBase, used only by env's audio-rate
// call site (gonnaDo). Identical math to cheapCurveBase; the per-shape
// sub-expressions are passed in precomputed instead of recomputed per sample:
//   sqrtInner    = sqrt(1/c - 1)
//   invSqrtInner = 1/sqrtInner
//   halfInvC     = 1/(2c)
// Drops 1 sqrt + 1 div (the 1/c) per sample bit-for-bit; the other two divisions
// become multiply-by-reciprocal, i.e. equal to cheapCurveBase to <=1 ULP (~1e-16,
// inaudible). For strict bit-exactness write atan(x/sqrtInner) and /(2*c) here
// instead (keeps the sqrt + 1/c win, leaves those two as divides).
cheapCurveBaseH(c, sqrtInner, invSqrtInner, halfInvC, x) = (log(c*(x*x-1)+1)+2*sqrtInner*atan(x*invSqrtInner)-2*x)*halfInvC;

cheapCurveReleaseS(c, invScale, zero, x) = (cheapCurveBase(c, x)-zero)*invScale;
cheapCurveAttackS(c, invScale, zero, x) = cheapCurveReleaseS(c, invScale, zero, x*-1+1)*-1+1;

// Base (unscaled) derivatives — these are the raw x(1-x)/denominator without
// dividing by curveScale.
derivativeBaseRelease(c, x) = x*(1-x)/(c*x*x+1-c);
derivativeBaseAttack(c, x) = x*(1-x)/(c*x*x+1-(2*c*x));

// Scaled derivatives (divided by curveScale) — used only where the actual
// envelope value is needed, not in the inverse path.
derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

// Variants taking precomputed 1/curveScale. Same result as derivativeRelease/Attack.
derivativeReleaseS(c, invScale, x) = derivativeBaseRelease(c, x)*invScale;
derivativeAttackS(c, invScale, x) = derivativeBaseAttack(c, x)*invScale;

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
// Phase of the release derivative's peak. The release base derivative is the
// horizontal flip of the attack one (derivativeBaseAttack(c,x) ==
// derivativeBaseRelease(c,1-x)), so its peak is the mirror of the attack peak.
peakPhaseRelease(c) = 1-peakPhaseAttack(c);

// Unscaled peak value — avoids dividing by curveScale then multiplying back
maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));
// Scaled version kept for reference / other uses
maxDerivativeAttack(c) = derivativeAttack(c, peakPhaseAttack(c));

inverseDerivativePart(c, D) = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));

inverseDerivativeTopRelease(c, D) = (1+inverseDerivativePart(c, D))/(2*(c*D+1));

inverseDerivativeBottomRelease(c, D) = (1-inverseDerivativePart(c, D))/(2*(c*D+1));

// Attack inverses via the same horizontal-flip relationship
inverseDerivativeTopAttack(c, D) = 1-inverseDerivativeTopRelease(c, D);

inverseDerivativeBottomAttack(c, D) = 1-inverseDerivativeBottomRelease(c, D);

// Original scaled versions kept for reference
// inverseDerivativePart(c, x) = sqrt(max(0, 1-4*(curveScale(c)*x-c*curveScale(c)*x)*(c*curveScale(c)*x+1)));
// inverseDerivativeTopRelease(c, x) = (1+inverseDerivativePart(c, x))/(2*(c*curveScale(c)*x+1));
// inverseDerivativeBottomRelease(c, x) = (1-inverseDerivativePart(c, x))/(2*(c*curveScale(c)*x+1));
// inverseDerivativeTopAttack(c, x) = inverseDerivativeTopRelease(c, x)*-1+1;
// inverseDerivativeBottomAttack(c, x) = inverseDerivativeBottomRelease(c, x)*-1+1;

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));

// ============================================================================
//  AUC level compensation — formula (SOURCE OF TRUTH) + superseded LUT
// ----------------------------------------------------------------------------
//  aucLevelMultFormula is the exact definition. AUC_attack(c) has no closed
//  form, so it is a midpoint-rule numerical integral (aucM steps). This is the
//  single source of truth: auc_poly.lib (imported at the top) is a polynomial
//  fit of it. Nothing in the audio path calls these definitions, so Faust drops
//  them from the compiled plugin — they exist only so regen_auc_poly.sh can
//  read aucLevelMultFormula via library() and re-bake the polynomial.
//
//  Everything below aucLevelMultFormula (aucTblContent, aucReadRaw,
//  table_aucLevelMult) is the older LUT path, kept for reference. It is dead
//  code and needs auc_waveform.lib (commented out at the top) re-enabled to
//  compile, since aucTbl lives there.
//
//  NOTE: compile with -double. At very small c (shape slider near 0) the
//  cheapCurveBase formula suffers catastrophic cancellation in single
//  precision, shifting the low-shape factor by ~0.8%. -double makes it exact.
// ============================================================================
aucN = 257;
// table points across the shape slider range [0,1]
aucM = 128;
// midpoint integration steps (init-time only; M=32 -> AUC err ~2e-7)

// AUC_attack reduced algebraically to a SINGLE integral of cheapCurveBase:
//   cheapCurveAttack(c,x) = 1 - (base(c,1-x)-base(c,0))/scale,  scale=base(c,1)-base(c,0)
//   => AUC_attack = 1 - (Ibase - base(c,0))/scale,  Ibase = integral_0^1 base(c,u) du
// This folds ~4x faster than integrating the fully composed curve (which would
// recompute scale and the endpoints at every sample point).
aucBaseIntegral(c) = sum(k, aucM, cheapCurveBase(c, (k+0.5)/aucM))/aucM;
aucAttack(c) = 1-(aucBaseIntegral(c)-cheapCurveBase(c, 0))/(cheapCurveBase(c, 1)-cheapCurveBase(c, 0));
aucMinAttack = aucAttack(shapeMap(1.0));
// smallest AUC -> normaliser
aucLevelMultFormula(s) = aucMinAttack/aucAttack(shapeMap(s));

// Two alternative compile paths, kept for reference — both compile too slowly:
//   1. aucReadRaw via rdtable(aucN, aucTblContent, idx) (commented out just
//      below): the integral as LIVE rdtable content makes Faust expand
//      aucN*aucM transcendentals at compile time.
//   2. the baked waveform (auc_waveform.lib): folding its read into att/rel
//      re-expands the 257-entry table at every att/rel use site (size *
//      fan-out).
// The live path uses the auc_poly.lib polynomial. Regenerate with:
//   regen_auc_poly.sh shapedSmoother.dsp 8 257
aucTblContent = aucLevelMultFormula(float(ba.time)/float(aucN-1)):min(1);
// aucReadRaw(idx) = rdtable(aucN, aucTblContent, idx);
aucReadRaw(idx) = aucTbl, idx:rdtable;
table_aucLevelMult(s) = it.interpolate_linear(frac, lo, hi):max(0):min(1)
    with {
        pos = s:max(0):min(1):*(aucN-1);
        i0 = int(pos);
        i1 = min(i0+1, aucN-1);
        frac = pos-i0;
        lo = aucReadRaw(i0);
        hi = aucReadRaw(i1);
    };

//============================================================================
// DJ-style "clip the peaks" gain stage  (mono, feedback topology)
//
// Two-state feedback loop produces, per sample:
//   gain : FAST gain-reduction envelope (instant attack, adaptive release)
//   ref  : SLOW reference envelope that trails `gain`
// clipPeaks then crossfades ref -> gain by how far the fast envelope has
// dipped below ref, and the result is scaled by `strength`.
//
// NB: all smoothing is in the LINEAR domain. A *linear* gain RISES when GR is
// RELEASED, so the att/rel time args to the smoothers are inverted vs the dB
// intuition. That is deliberate.
//============================================================================
DJcompression_gain_mono(strength, thresh, att, rel, knee, level) = loop~(_, _)// feedback states: prevGain, prevRef
:clipPeaks*strength
    with {
        // ref + (gain-ref)*clipMix  ==  lerp(ref -> gain).
        // clipMix rises 0->1 as the *current* fast GR sinks refBot -> refTop:
        // small peaks ride `ref`, big peaks follow the fast `gain`.
        clipPeaks(gain, ref) = ref+fastGRnow*clipMix
            with {
                fastGRnow = gain-ref;
                clipMix = it.remap(refBot, refTop, 0, 1, clamp(refBot, refTop, fastGRnow));
            };

        loop(prevGain, prevRef) = gain, ref
            with {
                fastGR = prevGain-prevRef;
                // fast GR the loop sees (prev sample)

                // FAST envelope: instant attack; release time set adaptively below.
                gain = gain_computer(1, thresh, knee, level):ba.db2linear:smootherARorder(maxOrder, orderRel, orderAtt, adaptiveRel, 0):ba.linear2db;

                // deeper fast GR -> longer release (fade_to_inf blows up as relWeight->1).
                adaptiveRel = fade_to_inf(relWeight, rel)*1@fullLatency;
                relWeight = 1-it.remap(dvBot, dvTop, 1, 0, clamp(dvBot, dvTop, fastGR));

                // SLOW reference envelope: 1st-order smoother trailing prevGain.
                ref = (prevGain-dvBot):ba.db2linear:smootherOrder(1, 1, refRel, 0):ba.linear2db;

                // ref release time: log-interpolated between slowRelease and ~infinity.
                refRel = interpolate_logarithmic(refWeight, slowRelease, slowRelease/ma.EPSILON)*1@fullLatency;
                refWeight = it.remap(refBot, refTop, 1, 0, clamp(refBot, refTop, fastGR));

                clamp(lo, hi, x) = x:min(hi):max(lo);
            };
    };
