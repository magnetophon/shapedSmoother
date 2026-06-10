declare name "shapedSmoother_core";
declare version "0.4";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — core algorithm, v0.4
//
//  CHANGES vs v0.3 (after field testing on the stepped test signal):
//
//  A. GUARDED LANDING. The incremental span (Step 3) trades the old loose
//     position-consistency of the span for noise immunity; under
//     mid-transition retargets the (phase, span, position) triple can
//     therefore drift apart, and the phase can complete with a real gap
//     left to the target. The v0.3 landing snapped across it regardless —
//     a visible teleport (up to ~7% of full scale measured on the stepped
//     test signal). The snap now only absorbs residuals up to one
//     peak-speed sample (anything the curve itself could have covered in a
//     sample); larger residuals reset the state WITHOUT snapping, so the
//     leftover distance gets a fresh shaped transition that velocity-
//     matches into the current speed instead of jumping or restarting from
//     zero.
//
//  B. CHEAP INCREMENTS. v0.3 computed the increment as an exact
//     antiderivative difference: two extra log+atan per sample over v0.1.
//     Bit-exactness of the increment only mattered while velocity matching
//     re-derived the phase from the realized speed every sample (v0.1);
//     with the trusted phase it just has to be good enough that the
//     landing absorbs the residual. Two-point Gauss-Legendre quadrature on
//     the cheap rational derivative has per-transition error ~1e-10 of the
//     span — five orders below the v0.1 one-point rule that caused the
//     original kink — at v0.1's CPU cost.
//
//  Merge of the two v0.2 lineages, keeping the best half of each:
//
//  FROM THE "exact" BRANCH (lands on target, on time):
//
//  1. EXACT INCREMENTS. The per-sample delta is the exact difference of the
//     curve's antiderivative: delta = totalStep*(fd(newPhase)-fd(anchor)).
//     Zero quadrature error; over an undisturbed transition the increments
//     telescope to exactly totalStep. (The v0.1 lineage re-derived the phase
//     from the measured speed every sample, which pins the integrator to
//     right-endpoint point samples of the derivative; the resulting drift
//     produced the end-of-transition kink, and even with the projection
//     corrected, transitions landed ~17% early.)
//
//  2. TRUSTED PHASE ADVANCE. The stored phase is the truth on samples where
//     nothing changed; velocity matching is demoted to an exception handler
//     for genuine discontinuities (target moved, direction flipped, outside
//     modification of the envelope — detected via the prevExpected state).
//
//  3. MIDPOINT-CORRECTED MATCHING. With exact increments the measured speed
//     is the MEAN derivative over the last step; the inverse-derivative
//     functions invert the POINT derivative and so return (approximately)
//     the midpoint of that step. Adding step/2 converts that to the current
//     position.
//
//  FROM THE noise-robustness BRANCH:
//
//  4. INCREMENTAL totalStep. The span from the (fixed) start of a transition
//     to the (moving) target changes only when the target moves, and by
//     exactly how much it moves: totalStep = prevTotalStep+(lookaheadX-
//     lookaheadX'). No needsRecompute heuristic, no round-trip of the state
//     through the curve, nothing for input noise to pump. This REPLACES the
//     exact branch's entire Step-5b monotone-anchor/effSpan machinery, which
//     existed to survive the per-sample span recompute under a moving
//     target; with the feedback removed at the source, that machinery is
//     dead weight (and measurably added spurious motion under noise).
//     It also fixes a real bug in the old entry path: rawTotalStep evaluated
//     fracDone with the NEW direction's curve on the OLD direction's stored
//     phase, so a direction flip with stale state could poison the new span
//     (e.g. a release born with the attack's span of -0.4).
//
//  5. SUPPORT FLOOR + LOWER OUTPUT CLAMP. The span may shrink only down to
//     the smallest span whose curve can still carry the current speed;
//     below it the envelope would slow discontinuously when clampedRatio
//     saturates. On the floor it decelerates along the curve at the curve's
//     steepest rate; the deliberate overshoot this allows is caught by the
//     lower output clamp (the floor and the clamp are a pair: the floor
//     without its catcher produced 0.066 of attack overshoot in testing).
//
//  6. FLOAT32 CONDITIONING. The curve formulas are rewritten as sums of
//     non-negative terms (see the note above cheapCurveBase); the naive
//     forms lose ~3.5 of single precision's ~7 digits at sharp shapes.
//
//  AND ONE FIX FOUND IN THE MERGE:
//
//  7. ROBUST LANDING. The landing test was newPhaseRaw >= 1, but the final
//     newPhaseRaw is (1-step)+step, which in floating point can round just
//     below 1 (it does for step = 1/480). The snap then never fires, the
//     one-epsilon overshoot of the exact increments flips the direction
//     with stale state, and the envelope glitches and creeps. The phase
//     advances in whole steps, so half a step of slack is unambiguous:
//     landed = active & (newPhaseRaw >= 1-halfStep).
// ============================================================================

import("stdfaust.lib");

// ============================================================================
//  GUI
// ============================================================================
MainGroup(x) = hgroup("[0]shapedSmoother", x);
TestGroup(x) = vgroup("[0]Test signal", x);
SmootherGroup(x) = vgroup("[1]Smoother", x);

// --- Test signal ---
testNoiseLevel = TestGroup(hslider("[0]noise level", 0, 0, 1, 0.001));
testNoiseRate = TestGroup(hslider("[1]noise rate", 42, 1, 1000, 1));
testBlockscale = TestGroup(hslider("[2]blockscale", 1, 0.01, 10, 0.01));
testFreq = TestGroup(hslider("[3]freq", 1, 0.001, 30, 0.001));
testStep1 = TestGroup(hslider("[4]step1", 0.75, -1, 1, 0.001));
testStep2 = TestGroup(hslider("[5]step2", 0.125, -1, 1, 0.001));
testSelect = TestGroup(checkbox("[6]signal select"));

testSignal = select2(testSelect, testSignal1, testSignal2);

testSignal1 = os.lf_squarewave(testFreq)*0.5;
testSignal2 = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };

// ============================================================================
//  THE SHAPED CURVES
// ============================================================================

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));

//  NOTE on arithmetic: the log argument is written c*x*x+(1-c) rather than
//  the algebraically identical c*(x*x-1)+1. The latter computes a tiny
//  result (~1-c at small x) as the difference of two near-unit numbers —
//  in single precision at sharp shapes (1-c ~ 3e-4) that cancellation
//  costs ~3.5 of float32's ~7 digits, and the noise lands on everything
//  built from the curve. The rewritten form is a sum of non-negative
//  terms: fully conditioned. (1-c) itself is exact in float for c >= 0.5
//  (Sterbenz), i.e. precisely the sharp shapes that need it.

cheapCurveBase(c, x) = (log(c*x*x+(1-c))+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);

curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x) = 1-cheapCurveRelease(c, 1-x);

//  Denominators rewritten as sums of non-negative terms for float32
//  conditioning at sharp shapes (see the note above cheapCurveBase):
//    c*x*x+1-c       == c*x*x+(1-c)
//    c*x*x+1-2*c*x   == c*(x-1)*(x-1)+(1-c)
derivativeBaseRelease(c, x) = x*(1-x)/(c*x*x+(1-c));
derivativeBaseAttack(c, x) = x*(1-x)/(c*(x-1)*(x-1)+(1-c));

derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
peakPhaseRelease(c) = 1-peakPhaseAttack(c);

maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));

inverseDerivativePart(c, D) = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));

inverseDerivativeTopRelease(c, D) = (1+inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeBottomRelease(c, D) = (1-inverseDerivativePart(c, D))/(2*(c*D+1));

inverseDerivativeTopAttack(c, D) = 1-inverseDerivativeTopRelease(c, D);
inverseDerivativeBottomAttack(c, D) = 1-inverseDerivativeBottomRelease(c, D);

// ============================================================================
//  THE ENVELOPE FOLLOWER
// ============================================================================

// --- Smoother controls (raw) ---
attMs = SmootherGroup(hslider("[0]att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01));
attackShapeSlider = SmootherGroup(hslider("[1]attack shape", 0, 0, 1, 0.001));
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));

// Derived parameters
maxHold = 0.05;
maxSR = 48000;
maxHoldSamples = maxHold*maxSR;

att = attMs/1000;
rel = relMs/1000;
att_samples = att*ma.SR:max(1);
rel_samples = rel*ma.SR:max(1);

process = testSignal<:(_@int(att_samples), shapedSmoother(_));

shapedSmoother(x) = lookaheadX, delayedX:env~(_, _, _, _):(_, _, _, !)
    with {
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHoldSamples);
        delayedX = x@int(att_samples);

        // --- Precomputed per-shape constants (slider-rate) ---
        attackShape = shapeMap(attackShapeSlider);
        releaseShape = shapeMap(releaseShapeSlider);

        attackInvCurveScale = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackZero = cheapCurveBase(attackShape, 0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
        releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

        // Fencepost: a transition occupies the CLOSED interval of N+1 output
        // samples [entry .. entry+N], so the phase advances in N+1 steps of
        // 1/(N+1). With 1/N the attack landed one sample early: lookaheadX
        // drops at m (the new minimum enters the window immediately) while
        // delayedX delivers the transient at m+480, and N increments
        // starting ON the entry sample complete at m+479. With 1/(N+1) the
        // attack's landing snap fires on exactly the sample where delayedX
        // presents the transient, and a release elapses exactly
        // rel_samples from onset to target.
        attStep = 1/((att*ma.SR):max(1)+1);
        relStep = 1/((rel*ma.SR):max(1)+1);

        // =======================================================================
        //  env: the recursive core
        //
        //  State carried between samples:
        //    prev          — current envelope value
        //    prevPhase     — where we are on the curve [0, 1]
        //    prevTotalStep — total span of the current transition
        //    prevExpected  — the value THIS loop computed last sample.
        //                    prev != prevExpected means something outside the
        //                    loop (the output clamps, or an external
        //                    algorithm) modified the envelope, so the stored
        //                    phase can no longer be trusted and we must
        //                    velocity-match.
        // =======================================================================
        env(prev, prevPhase, prevTotalStep, prevExpected, lookaheadX, delayedX) = result, newPhase, totalStepOut, expected
            with {
                // --- Step 1: Are we attacking, releasing, or idle? ---
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                // --- Step 2: Select shape parameters for current direction ---
                shape = select2(releasing, attackShape, releaseShape);
                invCurveScale = attackInvCurveScale+releasing*(releaseInvCurveScale-attackInvCurveScale);
                zeroVal = select2(releasing, attackZero, releaseZero);
                maxDerivVal = select2(releasing, attackMaxDerivBase, releaseMaxDerivBase);
                step = select2(releasing, attStep, relStep);
                halfStep = 0.5*step;

                // --- Curve helpers ---
                // cbAt(p): normalized curve value (attack flip applied inside)
                // fdOf(p): fraction of the transition completed at phase p,
                //          0 at the start, 1 at the target, for BOTH directions.
                //          Position on curve = transitionStart + totalStep*fdOf(p).
                cbAt(p) = (cheapCurveBase(shape,
                    select2(releasing, 1-p, p))-zeroVal)*invCurveScale;
                fdOf(p) = select2(releasing, 1-cbAt(p), cbAt(p));

                // --- Step 3: totalStep (incremental) ---
                //
                // totalStep is the distance from the (fixed) start of the
                // current transition to the (moving) target. The envelope's
                // own progress lives in the phase, not here — so the only
                // thing that can change the span mid-transition is the
                // target moving. Integrate exactly that:
                //
                //   totalStep = prevTotalStep + (lookaheadX - lookaheadX')
                //
                // A static target gives a zero increment (an exact hold); a
                // moving target shrinks or grows the span by exactly its own
                // movement, in either direction. The state never round-trips
                // through the curve, so there is no accounting drift to
                // re-inject and nothing for input noise to pump: zero-mean
                // target jitter sums to the target's net movement by
                // construction.
                //
                // On entry into a direction the span is the full remaining
                // distance. Entry is detected by the sign of prevTotalStep
                // (positive release, negative attack, zero idle) — never by
                // evaluating the curve on the previous direction's phase.
                //
                // supportTotalStep is the smallest span whose curve can
                // carry the current speed (the curve derivative maxes out at
                // maxDerivVal). When the target jumps closer mid-transition
                // the span may shrink only down to this floor; below it,
                // clampedRatio would saturate and the envelope would slow
                // discontinuously. On the floor it instead decelerates along
                // the curve at the curve's steepest rate, and the deliberate
                // overshoot this allows is caught by the output clamps in
                // Step 7. One expression serves both directions: it floors
                // the release span (both positive) and ceils the attack span
                // (both negative).
                incrementalTotalStep = prevTotalStep+(lookaheadX-lookaheadX');
                entryTotalStep = lookaheadX-prev;
                supportTotalStep = prevSpeed/(maxDerivVal*invCurveScale*step);

                totalStep = select2(releasing,
                    min(select2(prevTotalStep>=0, incrementalTotalStep, entryTotalStep),
                        supportTotalStep),
                    max(select2(prevTotalStep<=0, incrementalTotalStep, entryTotalStep),
                        supportTotalStep))*active;

                // --- Step 4: trusted phase advance vs. velocity matching ---
                //
                // On a sample where nothing changed, the stored phase IS the
                // truth: re-deriving it from the measured speed can only add
                // error. So velocity matching is used only when:
                //   - the target moved,                  (sameTarget fails)
                //   - the direction flipped,             (sameDirection fails)
                //   - the envelope was modified outside  (undisturbed fails)
                //     this loop (external algorithms, or our own output
                //     clamps having altered the previous output),
                //   - or we are entering a transition.   (prevPhase == 0)
                sameTarget = (lookaheadX==lookaheadX');
                sameDirection = (releasing==releasing')&(attacking==attacking');
                undisturbed = (prev==prevExpected);
                trusted = sameTarget&sameDirection&undisturbed&(prevPhase>0);

                // --- Step 5: velocity matching (for untrusted samples) ---
                //
                // Since the emitted delta is an exact curve increment, the
                // measured speed is the MEAN derivative over the last step.
                // The inverse-derivative functions invert the POINT
                // derivative, so they return (approximately) the midpoint of
                // that step. Adding halfStep converts the recovered midpoint
                // into the current position. (Without this, the phase clock
                // gains ~step/2 per matched sample.)
                //
                // The correction only makes sense when a previous interval
                // in the same direction actually exists; on a fresh entry
                // into a transition (speed ~0, no history) it would skip the
                // first half-step of curve and leave a small residual to be
                // absorbed by the landing snap.
                matchCorr = halfStep*(sameDirection&(prevPhase>0));
                prevSpeed = prev-prev';
                speedRatio = prevSpeed/(totalStep*invCurveScale*step+(1-active)*1e-30);
                clampedRatio = max(0, min(speedRatio, maxDerivVal));

                gonnaDo(phase) = (1-fdOf(phase))*totalStep;

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio))+matchCorr:min(1))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                matchedPos = select2(releasing,
                    // Attack branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    // Release branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        inverseDerivativeTopRelease(shape, clampedRatio)))+matchCorr:max(0):min(1-step);

                // (The old Step 5b — the monotone anchor, spanContinuous and
                // effSpan machinery — is gone. It compensated for the
                // per-sample span recompute feeding back into the speed
                // ratio under a moving target; the incremental totalStep of
                // Step 3 removes that feedback at the source, and with it
                // the machinery measurably ADDED spurious motion under
                // noise.)

                // --- Step 6: advance phase ---
                //
                // anchor   = where we are on the curve NOW
                // newPhase = where we will be after this sample
                //
                // The 1-step clamp on the anchor makes the phase reach 1 (to
                // within rounding) on the final sample of a transition; the
                // landing test in Step 7 has half a step of slack to absorb
                // that rounding.
                anchor = select2(trusted, matchedPos, prevPhase:min(1-step));
                newPhaseRaw = (anchor+step)*active;

                // --- Step 7: curve increment + guarded landing ---
                //
                // The increment integrates the curve derivative over
                // [anchor, anchor+step] with 2-point Gauss-Legendre
                // quadrature on the cheap rational derivative (see header
                // note B; the derivative of fdOf is the same expression for
                // both directions — the attack flip's two sign changes
                // cancel in the chain rule).
                glA = 0.21132486540518713;
                glB = 0.78867513459481287;
                derivOf(p) = derivativeBaseRelease(shape,
                    select2(releasing, 1-p, p))*invCurveScale;
                delta = totalStep*step*0.5*(derivOf(anchor+glA*step)+derivOf(anchor+glB*step));

                // The value this loop wants to write. Stored as state so the
                // next sample can detect outside modification (see Step 4).
                expected = prev+delta;

                // Landing: when the phase completes, snap to the target —
                // but only across a residual the snap was designed for
                // (quadrature epsilon, float epsilon). guardSize is one
                // peak-speed sample. Under mid-transition retargets the
                // phase can complete with a real gap left (see header note
                // A); snapping across that is a teleport, so instead the
                // state resets WITHOUT snapping and the next sample starts
                // a fresh shaped transition over the leftover distance,
                // entered through velocity matching so it continues from
                // the current speed.
                //
                // The phase-completion test needs half a step of slack: the
                // final newPhaseRaw is (1-step)+step, which can round just
                // below 1 (it does for step = 1/480), and a missed landing
                // flips the direction on the residual epsilon with stale
                // state. The phase advances in whole steps, so halfStep
                // slack is unambiguous.
                //
                // The output never goes above the raw (delayed) input and
                // never below the attack target. The lower clamp is the
                // required catcher for supportTotalStep's deliberate
                // overshoot (Step 3); it is inert outside of attack, since
                // min(prev, lookaheadX) equals prev while releasing or
                // idle, and it cannot conflict with the upper clamp because
                // lookaheadX <= delayedX by construction (the sliding-min
                // window contains the delayed sample).
                phaseDone = active&(newPhaseRaw>=(1-halfStep));
                guardSize = abs(totalStep)*maxDerivVal*invCurveScale*step;
                landed = phaseDone&(abs(lookaheadX-expected)<=guardSize);

                result = select2(landed,
                    min(expected, delayedX),
                    min(lookaheadX, delayedX)):max(min(prev, lookaheadX));

                // Completion (landed or not) fully RESETS the transition
                // state. With the state cleared: constant target -> idle as
                // before; drifting target -> fresh micro-transitions that
                // follow at curve speed; big jump -> a full new
                // attack/release; un-snapped residual -> a fresh shaped
                // transition over the leftover.
                newPhase = newPhaseRaw*(1-phaseDone);
                totalStepOut = totalStep*(1-phaseDone);
            };
    };

// ============================================================================
//  GLOSSARY
// ============================================================================
//
//  phase:
//    A value in [0, 1] tracking progress through the current attack or
//    release transition. 0 = just started, 1 = done. Advances by `step`
//    per sample; trusted between discontinuities, re-derived from the
//    measured speed (velocity matching) across them.
//
//  totalStep:
//    The full amplitude span (in linear gain) of the current transition,
//    from where it started to where it's heading. Maintained incrementally:
//    held exactly while the target is static, moved by exactly
//    (lookaheadX - lookaheadX') when it isn't, recomputed as
//    (lookaheadX - prev) on entry into a direction, and never allowed to
//    shrink past the span whose curve can still carry the current speed
//    (supportTotalStep).
//
//  step:
//    Phase increment per sample = 1 / (duration_in_samples).
//
//  fdOf(p):
//    Fraction of the transition completed at phase p, for both directions.
//    Position on the curve = transitionStart + totalStep*fdOf(p). The
//    per-sample delta integrates its derivative over one step with 2-point
//    Gauss-Legendre quadrature; the ~1e-10 per-transition residual is
//    absorbed by the landing snap.
//
//  velocity matching:
//    Across a discontinuity (retarget, direction flip, outside
//    modification), find the phase on the (possibly new) curve whose
//    derivative equals the envelope's current speed, then continue from
//    there. With exact increments the measured speed is a mean over the
//    last step, hence the halfStep matchCorr.
//
//  gonnaMakeIt:
//    "If we re-anchor on the decelerating branch now, will the resulting
//    trajectory still reach the target?" Picks the accelerating vs
//    decelerating branch of the inverse derivative.
//
//  landing:
//    When the phase completes, the transition state is fully reset, and
//    the output snaps to the target only if the residual is at most one
//    peak-speed sample (quadrature/float epsilon territory). Larger
//    residuals — the (phase, span, position) drift that mid-transition
//    retargets cause under the incremental span — are NOT snapped
//    (teleport); the leftover gets a fresh shaped transition that
//    velocity-matches into the current speed. The phase-completion test
//    carries halfStep of slack because the final phase value is produced
//    by accumulation and can round below 1.
//
//  lookaheadX:
//    The sliding minimum of the input over the next att_samples: the
//    envelope's target.
//
//  AUC compensation:
//    (Omitted from this core file for clarity.) Adjusts the attack/release
//    duration so that changing the shape doesn't change the average level
//    of the envelope — sharper shapes get proportionally shorter durations.
