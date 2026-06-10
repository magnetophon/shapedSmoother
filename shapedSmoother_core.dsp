declare name "shapedSmoother_core";
declare version "0.2";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — core algorithm, v0.2
//
//  CHANGE vs v0.1 (the "end-of-release kink" fix):
//
//  The old loop integrated the curve with a one-point rectangle rule
//  (delta = derivative(newPhase)*step). On the decelerating side this
//  systematically undershoots, so the envelope's amplitude drifted behind
//  the phase. Near the end of a sharp release, gonnaMakeIt would flip to
//  "won't make it", hop to the accelerating branch, and sprint the residual
//  gap: a visible kink, and an effective release time much shorter than
//  `rel`.
//
//  The fix has three parts:
//
//  1. EXACT INCREMENTS. The per-sample delta is now the exact difference of
//     the curve's antiderivative: delta = totalStep*(fd(newPhase)-fd(anchor)).
//     Position can no longer drift off the curve, so the end-of-transition
//     catch-up never fires.
//
//  2. TRUSTED PHASE ADVANCE. The phase is only re-derived from the measured
//     velocity when something actually changed (target moved, direction
//     flipped, or `prev` was modified by an outside clamp/algorithm,
//     detected by comparing prev against the value this loop wrote last
//     sample). On undisturbed samples the phase simply advances by `step`
//     from the stored prevPhase, which keeps both amplitude AND timing
//     exact. Velocity matching is reserved for genuine discontinuities,
//     where its small anchoring error is masked by the transition itself.
//
//  3. MIDPOINT-CORRECTED MATCHING. When velocity matching IS used: with
//     exact increments the measured speed is the MEAN derivative over the
//     last step, not the point derivative, so the inverse-derivative
//     functions return (approximately) the midpoint of that step.
//     Adding step/2 converts that to the current position. Without this,
//     the phase clock gains ~half a step per matched sample.
//
//  A transition now also LANDS exactly: when the phase reaches 1, the
//  output snaps to the target (this absorbs the last floating-point
//  epsilon instead of leaving the loop creeping forever).
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
testFreq = TestGroup(hslider("[3]freq", 1, 0.001, 4, 0.001));
testStep1 = TestGroup(hslider("[4]step1", 0.75, -1, 1, 0.001));
testStep2 = TestGroup(hslider("[5]step2", 0.125, -1, 1, 0.001));

testSignal = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };

// ============================================================================
//  THE SHAPED CURVES (unchanged)
// ============================================================================

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));

cheapCurveBase(c, x) = (log(c*(x*x-1)+1)+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);

curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x) = 1-cheapCurveRelease(c, 1-x);

derivativeBaseRelease(c, x) = x*(1-x)/(c*x*x+1-c);
derivativeBaseAttack(c, x) = x*(1-x)/(c*x*x+1-(2*c*x));

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

        releasePeakPhase = peakPhaseRelease(releaseShape);

        attStep = 1/((att*ma.SR):max(1));
        relStep = 1/((rel*ma.SR):max(1));

        // =======================================================================
        //  env: the recursive core
        //
        //  State carried between samples:
        //    prev          — current envelope value
        //    prevPhase     — where we are on the curve [0, 1]
        //    prevTotalStep — total span of the current transition
        //    prevExpected  — the value THIS loop computed last sample.
        //                    prev != prevExpected means something outside the
        //                    loop (the delayedX ceiling, or an external
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

                cbPrev = cbAt(prevPhase);
                fracDone = select2(releasing, 1-cbPrev, cbPrev);

                // --- Step 3: totalStep (unchanged from v0.1) ---
                rawTotalStep = (lookaheadX-prev)+prevTotalStep*fracDone;
                needsRecompute = (lookaheadX!=lookaheadX')|(prevTotalStep<=0);

                totalStep = select2(releasing,
                    rawTotalStep:min(prevTotalStep),
                    select2(needsRecompute, prevTotalStep, rawTotalStep))*active;

                // --- Step 4: trusted phase advance vs. velocity matching ---
                //
                // On a sample where nothing changed, the stored phase IS the
                // truth: re-deriving it from the measured speed can only add
                // error. So velocity matching is now used only when:
                //   - the target moved,                  (sameTarget fails)
                //   - the direction flipped,             (sameDirection fails)
                //   - the envelope was modified outside  (undisturbed fails)
                //     this loop (external algorithms, or our own delayedX
                //     ceiling having clamped the previous output),
                //   - or we are entering a transition.   (prevPhase == 0)
                sameTarget = (lookaheadX==lookaheadX');
                sameDirection = (releasing==releasing')&(attacking==attacking');
                undisturbed = (prev==prevExpected);
                trusted = sameTarget&sameDirection&undisturbed&(prevPhase>0);

                // --- Step 5: velocity matching (for untrusted samples) ---
                //
                // Same idea as before, with one correction: since the emitted
                // delta is now an exact curve increment, the measured speed is
                // the MEAN derivative over the last step. The inverse-
                // derivative functions invert the POINT derivative, so they
                // return (approximately) the midpoint of that step. Adding
                // halfStep converts the recovered midpoint into the current
                // position. (Without this, the phase clock gains ~step/2 per
                // matched sample.)
                //
                // The correction only makes sense when a previous interval in
                // the same direction actually exists; on a fresh entry into a
                // transition (speed ~0, no history) it would skip the first
                // half-step of curve and leave a small residual to be
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

                // --- Step 5b: monotone anchor during a continuing transition ---
                //
                // With a continuously moving target, the span recomputes
                // every sample and feeds back into the speed ratio:
                //   - span shrinking (closing target): ratio saturates above
                //     maxDerivVal, the clamp maps to the derivative's peak,
                //     and the anchor gets re-pinned there every sample;
                //   - span growing (receding target): the same speed maps to
                //     an ever-earlier phase, so the anchor is pushed
                //     backward and the envelope's acceleration is endlessly
                //     restarted, while the bottom branch caps the phase just
                //     below the peak.
                // Both show up as the phase buzzing at/under the peak while
                // the envelope tracks badly. The cure: while the transition
                // CONTINUES (same direction, same gonnaMakeIt branch), the
                // anchor may never move backward — and during pursuit
                // (gonnaMakeIt false) it parks AT the peak: maximum-speed
                // chase, which is the best the curve can do. Backward jumps
                // remain allowed exactly when they are legitimate: on a
                // branch change (e.g. the target leaps away near landing and
                // we must restart accelerating).
                // The rule must bind ONLY on recurring matching (the noise
                // feedback), never on the FIRST matched sample after a
                // trusted stretch: that first sample is the legitimate
                // velocity match of a discrete retarget, and constraining it
                // trades speed continuity for a visible kink.
                attackPeakPhase = peakPhaseAttack(attackShape);
                peakPhaseDir = select2(releasing, attackPeakPhase, releasePeakPhase);
                wasMatching = 1-trusted';
                monotone = wasMatching&sameDirection&(prevPhase>0)&(gonnaMakeIt==gonnaMakeIt');
                anchorMono = max(matchedPos, prevPhase:min(1-step)):min(select2(gonnaMakeIt, peakPhaseDir, 1-step));
                anchorMatched = select2(monotone, matchedPos, anchorMono);

                // --- Step 6: advance phase ---
                //
                // anchor   = where we are on the curve NOW
                // newPhase = where we will be after this sample
                //
                // The 1-step clamp on the anchor guarantees newPhase reaches
                // exactly 1 on the final sample of a transition (landing),
                // instead of asymptotically creeping below the target.
                anchor = select2(trusted, anchorMatched, prevPhase:min(1-step));
                newPhaseRaw = (anchor+step)*active;

                // --- Step 7: exact curve increment + landing ---
                //
                // delta is the exact difference of the antiderivative between
                // anchor and newPhase: integrating the curve with zero
                // quadrature error. Summed over a whole undisturbed
                // transition this telescopes to exactly totalStep, so the
                // envelope arrives on target, on time, with no catch-up kink.
                // On monotone samples the span is being recomputed every
                // sample (moving target), so increments of the stored span no
                // longer telescope to the true remaining distance — marching
                // the phase to 1 would end in a teleporting landing snap.
                // Instead, cover the corresponding fraction of the ACTUAL
                // remaining distance: each sample moves
                //   (X-prev) * (fd(new)-fd(anchor)) / (1-fd(anchor)),
                // which lands exactly on the target at phase 1, continuously.
                // (This form cannot blow up: the increment never exceeds the
                // remaining distance.)
                fdA = fdOf(anchor);
                fdN = fdOf(newPhaseRaw);
                deltaCurve = totalStep*(fdN-fdA);
                deltaMono = (lookaheadX-prev)*(fdN-fdA)/max(1e-30, 1-fdA);
                delta = select2(monotone&(1-trusted), deltaCurve, deltaMono);

                // The value this loop wants to write. Stored as state so the
                // next sample can detect outside modification (see Step 4).
                expected = prev+delta;

                // When the phase reaches 1 the transition is complete: land
                // exactly on the target, absorbing any last floating-point
                // epsilon. Otherwise: move along the curve, but never above
                // the raw (delayed) input.
                landed = active&(newPhaseRaw>=1);
                result = select2(landed,
                    min(expected, delayedX),
                    min(lookaheadX, delayedX));

                // Landing must fully RESET the transition state. If the
                // phase stayed parked at 1, an envelope tracking a slowly
                // drifting target through repeated landings would keep
                // monotone's anchor floor at 1 — and then a big target jump
                // would be swallowed by the landing snap (a teleport)
                // instead of restarting a proper shaped transition. With the
                // state cleared: constant target -> idle as before; drifting
                // target -> fresh micro-transitions that follow at curve
                // speed; big jump -> a full new attack/release.
                newPhase = newPhaseRaw*(1-landed);
                totalStepOut = totalStep*(1-landed);
            };
    };
