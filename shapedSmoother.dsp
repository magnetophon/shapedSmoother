declare name "shapedSmoother";
declare version "0.3.4";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

// v0.3.4 The LUT scales att/rel by aucLevelMult(shape) so changing the shape no longer
// changes the average level. See the "AUC level compensation" comment on the
// att/rel sliders and the "AUC level-compensation LUT" block lower down.
// Folding the factor into att/rel makes the reported latency shape-dependent
// (largest at the sharpest shape).
//
// v0.3.3: fix the release phase failing to reach 1 for some rel/release-shape
// combinations. The release branch chose its accelerate/decelerate phase via
// `gonnaMakeIt = projected > lookaheadX`, but for a single step `projected`
// equals the target by construction, so the test was decided by rounding. A
// spurious False past the derivative peak teleported the phase to the bottom
// branch (~0), permanently losing phase and making the output reach the target
// early (release faster than dialed). Replaced the release selector with a
// direct test — accelerate when before the derivative peak or when the target
// just grew (catch-up) — which matches gonnaMakeIt's intended value without the
// tie. Attack is unchanged (its analogous flip pins the phase near 1, benign)
// and the phase resync is unchanged.
//
// v0.3.2: v0.3.1 with the dead `remainingCurve` / `new_totalStep` code
// removed. In v0.3 Faust DCE'd it; in v0.3.1 (calling the new `*S` helpers)
// it was no longer DCE'd, adding ~5 transcendentals per sample to the inner
// loop. Removing the source-level dead code eliminates the regression.
//
// v0.3.1: v0.3 verbatim + slider-rate hoisting of per-shape constants.
// NO algebraic restructuring vs v0.3 — attack/release stay as separate
// function calls, derivativeBaseAttack/Release stay as separate functions,
// inverseDerivative*Attack/Release stay as 4 separate things. The only
// change is that quantities that depend solely on the shape sliders
// (which Faust failed to auto-lift through `select2(releasing, ...)`) are
// computed per-shape and selected scalar-only inside env().
//
// Specifically hoisted:
//   attackShape, releaseShape          (shapeMap of the slider)
//   attackInvCurveScale, releaseInvCurveScale  (1/curveScale)
//   attackZero, releaseZero            (cheapCurveBase(c, 0))
//   attackMaxDerivBase, releaseMaxDerivBase    (maxDerivativeBaseAttack)
//
// Everything else is byte-identical to v0.3.

process = test2@(att_samples+brake_samples+rel_hold_samples), shapedSmoother(test2);
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
        // These are pure functions of the shape sliders. Faust's auto-lift
        // doesn't push through `select2(releasing, sliderA, sliderB)` so in
        // v0.3 every `cheapCurveBase`/`curveScale`/`maxDerivativeBaseAttack`
        // call on `shape` fires audio-rate. Precompute them per-shape.
        attackShape = shapeMap(attackShapeSlider);
        releaseShape = shapeMap(releaseShapeSlider);

        attackInvCurveScale = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackZero = cheapCurveBase(attackShape, 0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
        releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

        // v0.3.3: phase of the release derivative's peak (slider-rate). Used to
        // pick the accelerating (bottom) vs decelerating (top) branch robustly.
        releasePeakPhase = peakPhaseRelease(releaseShape);

        env(prev, prevPhase, prevTotalStep, lookaheadX) = result, newPhase, totalStep
            with {
                // Audio-rate scalars: select the precomputed value, no math.
                shape = select2(releasing, attackShape, releaseShape);
                invCurveScale = select2(releasing, attackInvCurveScale, releaseInvCurveScale);
                zeroVal = select2(releasing, attackZero, releaseZero);
                maxDerivVal = select2(releasing, attackMaxDerivBase, releaseMaxDerivBase);

                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                activeTime = select2(releasing, att, rel);
                totalNRSteps = (activeTime*ma.SR):max(1);
                step = 1/totalNRSteps;

                // make delta dependent on how much GR we have to do until the full lookahead, so it self corrects

                todo = lookaheadX-prev;

                totalStep = select2(releasing,
                    (lookaheadX-prev):min(prevTotalStep),
                    (lookaheadX-prev):max(prevTotalStep))*active;

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
                        cb = (cheapCurveBase(shape, select2(releasing, 1-phase, phase))-zeroVal)*invCurveScale;
                    };

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio)))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                // v0.3.3: robust branch selection for the RELEASE.
                // For a single step, `projected` equals the target by
                // construction, so `gonnaMakeIt` is a floating-point coin-flip;
                // when it lands False past the derivative peak it sends the
                // phase to the BOTTOM branch (~0 for the release), losing phase
                // permanently and making the output race to target early.
                // The bottom branch is genuinely needed only (a) before the
                // derivative peak (the accelerating part, incl. release onset),
                // or (b) when the target has just grown and we must re-accelerate
                // to catch it. Detect both directly instead of via `gonnaMakeIt`:
                //   prevPhase < releasePeakPhase  -> still accelerating
                //   totalStep  > prevTotalStep    -> a larger target appeared
                // This reproduces gonnaMakeIt's *intended* value for a single
                // step without the tie, and preserves the mid-fade catch-up.
                // The phase resync (`phaseAtMatchingSpeed + step`) is unchanged.
                // ATTACK is left byte-identical: there the spurious flip pins
                // the phase near 1, where it is already heading, so it is benign.
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
                brakeMult = 1-brakeRamp;
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

// Raw shape sliders. Named here so the AUC level-compensation LUT (defined
// lower in the file) can index them, and so env() can shapeMap() them.
attackShapeSlider = hslider("attack shape", 0, 0, 1, 0.001);
releaseShapeSlider = hslider("release shape", 0, 0, 1, 0.001);

// AUC level compensation.
// Goal: changing the shape must NOT change the average level. A single-step
// response of duration T has time-area  T * AUC_attack(shape) * step. We scale
// the duration by  aucLevelMult(shape) = AUC_attack(shapeMap(1)) / AUC_attack(shape),
// so  T_eff * AUC = const  -> the area is shape-independent. The factor is
// normalized to <=1 (==1 at the sharpest shape), so durations only ever shorten.
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
// When off the multiplier is exactly 1.0, i.e. att/rel revert to the raw
// slider values — and the reported latency goes back to being
// shape-independent (see the AUC comment on the att/rel sliders).
attAucComp = checkbox("att auc comp");
relAucComp = checkbox("rel auc comp");

// Blend the multiplier toward 1.0 when the switch is off. Branchless on
// purpose: `on` is 0 or 1 so it's exact, the multiplier is clamped to [0,1]
// so `on*(m-1)` can't blow up, and this is slider-rate regardless.
aucLevelMultSwitched(on, s) = 1+on*(aucLevelMult(s)-1);

att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01)/1000*aucLevelMultSwitched(attAucComp, attackShapeSlider);
att_samples = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
rel = hslider("rel[scale:log]", 0.05*1000, 1, 5000, 0.1)/1000*aucLevelMultSwitched(relAucComp, releaseShapeSlider);
rel_samples = rel*ma.SR:max(1);

// Release hold (ported from lamb-rs): linear scale, allows 0 (off).
// Default 50 ms matches lamb's `[08]release hold`.
rel_hold = hslider("rel hold[unit:ms]", 50, 0, maxRelHold*1000, 1)*0.001;
rel_hold_samples = rel_hold*ma.SR;

// Brake window length per spec:
//   = 50 ms by default,
//   = max(att, rel) if both att AND rel are shorter than 50 ms.
// (min(50ms, max(att,rel)) gives exactly this.)
brake_samples = min(maxBrakeSamples, max(att_samples, rel_samples)):max(1);
test3 = os.lf_squarewave(hslider("freq", 1, 0.001, 4, 0.001));

// Test signal
test2 = it.interpolate_linear(hslider("noise level", 0, 0, 1, 0.001),
    (loop~_),
    no.lfnoise(hslider("noise rate", 42, 1, 1000, 1)))
    with {
        loop(prev) = no.lfnoise0(blockscale*(abs(prev*69)%9:pow(0.75)*5+1));
        blockscale = hslider("blockscale", 1, 0.01, 10, 0.01);
    };
testSig = os.lf_sawpos(0.5)<:(((_>0.25)*hslider("step1", 0.75, -1, 1, 0.001))+((_>0.5)*hslider("step2", 0.125, -1, 1, 0.001)));

// *************************************** the NEW curves: ******************************
// New curves by nuchi:
// https://www.wolframalpha.com/input?i=antiderivative+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29
// https://www.wolframalpha.com/input?i=inverse+function+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29+%2F+C
// https://www.desmos.com/calculator/9wtfhymvr0
// with scaling for c:
// https://www.desmos.com/calculator/dynyjjkuli
// with flip horizontal:
// https://www.desmos.com/calculator/bsr8cdn21v

cheapCurveBase(c, x) = (log(c*(pow(x, 2)-1)+1)+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);
curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);
cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x) = cheapCurveRelease(c, x*-1+1)*-1+1;

// v0.3.1: variants taking precomputed 1/curveScale and cheapCurveBase(c,0).
// Mathematically identical to cheapCurveRelease/Attack — they just don't
// recompute curveScale and cheapCurveBase(c,0) per sample.
//   invScale = 1/curveScale(c)
//   zero     = cheapCurveBase(c, 0)
cheapCurveReleaseS(c, invScale, zero, x) = (cheapCurveBase(c, x)-zero)*invScale;
cheapCurveAttackS(c, invScale, zero, x) = cheapCurveReleaseS(c, invScale, zero, x*-1+1)*-1+1;

// Base (unscaled) derivatives — these are the raw x(1-x)/denominator without
// dividing by curveScale.
derivativeBaseRelease(c, x) = x*(1-x)/(c*pow(x, 2)+1-c);
derivativeBaseAttack(c, x) = x*(1-x)/(c*pow(x, 2)+1-(2*c*x));

// Scaled derivatives (divided by curveScale) — used only where the actual
// envelope value is needed, not in the inverse path.
derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

// v0.3.1: variants taking precomputed 1/curveScale. Same math.
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
//  AUC level-compensation LUT
// ----------------------------------------------------------------------------
//  Reuses cheapCurveBase / shapeMap defined above. AUC_attack(c) has no
//  closed form, so it is a midpoint-rule numerical integral; that sum runs
//  only while the rdtable is filled at init (aucN*aucM evals once), never per
//  audio sample. Read with linear interpolation.
//
//  NOTE: compile with -double. At very small c (shape slider near 0) the
//  cheapCurveBase formula suffers catastrophic cancellation in single
//  precision, shifting the low-shape factor by ~0.8%. -double makes it exact.
// ============================================================================
aucN = 257;
// table points across the shape slider range [0,1]
aucM = 32;
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

aucTblContent = aucLevelMultFormula(float(ba.time)/float(aucN-1));
aucReadRaw(idx) = rdtable(aucN, aucTblContent, idx);

aucLevelMult(s) = it.interpolate_linear(frac, lo, hi):max(0):min(1)
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
