declare name "shapedSmoother_core";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — core algorithm
//  Extracted from shapedSmoother.dsp v0.3.6 by Bart Brouns (AGPL-3.0-only)
//
//  This file contains only the essential pieces:
//    1. The test signal generator
//    2. The shaped curve functions (attack & release shapes)
//    3. The envelope follower core (env)
//  with explanations of how everything fits together.
// ============================================================================

import("stdfaust.lib");

// ============================================================================
//  OVERVIEW
// ============================================================================
//
//  shapedSmoother is an envelope follower whose attack and release transitions
//  follow user-controllable curves rather than the usual exponential decay.
//
//  Signal flow (simplified):
//
//    input
//      |
//      v
//    slidingMin    — lookahead: finds the minimum in the next `att` samples,
//      |              so the envelope can START falling before a transient hits
//      v
//    env feedback  — the core: a recursive loop that moves `prev` (the current
//      |              envelope value) toward `lookaheadX` (the target) using a
//      |              shaped curve to determine the speed at each sample
//      v
//    output
//
//  The key insight is that attack and release are NOT exponential filters.
//  Instead, the envelope traverses a pre-defined curve shape over a fixed
//  number of samples. The curve shape is controlled by a single parameter
//  `c` (0..1) that morphs from s-shaped to very sharp/sudden.

// ============================================================================
//  1. TEST SIGNAL
// ============================================================================
//
//  A simple test signal for auditioning the smoother without external audio.
//  `lfnoise0` produces random steps; the output is a blend between that
//  deterministic chaos and a separate low-frequency noise source.
//
//  - testNoiseLevel: crossfade between the deterministic part and lfnoise
//  - testNoiseRate:  rate of the lfnoise component
//  - testBlockscale: scales the rate of the deterministic random-walk

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

testSignal = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        // Self-modulating random walk: the previous output determines
        // the rate of the next random step, giving organic variation.
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };

// ============================================================================
//  2. THE SHAPED CURVES
// ============================================================================
//
//  The core idea: instead of an exponential smoother, transitions follow
//  a curve defined by cheapCurveBase(c, x), where:
//    c  = shape parameter in (0, 1)   (0 = linear, 1 = extremely sharp)
//    x  = phase in [0, 1]             (progress through the transition)
//
//  The curve was designed by nuchi. It is the antiderivative of:
//      x * (1-x) / (c * x² + 1 - c)
//  which gives a smooth S-ish shape whose "sharpness" is controlled by c.
//
//  References:
//    https://www.desmos.com/calculator/bsr8cdn21v
//
//  IMPORTANT: the shape parameter `c` must be in (0, 1) exclusive — never
//  exactly 0 or 1, or we hit division by zero / log(0). shapeMap() below
//  maps the user-facing slider [0, 1] to a safe sub-range.

// --- 2a. shapeMap: slider → internal shape parameter ---
//
//  Maps the user's linear [0, 1] slider to an exponential curve that stays
//  safely inside (0, 1). At slider=0, shape ≈ 0.0001 (nearly linear);
//  at slider=1, shape ≈ 0.9999 (very sharp, nearly a step function).

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));

// --- 2b. cheapCurveBase: the raw (unscaled) curve ---
//
//  This is the antiderivative of the rational kernel x(1-x)/(cx²+1-c).
//  It contains a log term (the rational part) and an atan term (from the
//  partial-fraction decomposition). The result is NOT normalized to [0,1]
//  — that's what curveScale and the Release/Attack wrappers do.

//  NOTE on arithmetic: the log argument is written c*x*x+(1-c) rather than
//  the algebraically identical c*(x*x-1)+1. The latter computes a tiny
//  result (~1-c at small x) as the difference of two near-unit numbers —
//  in single precision at sharp shapes (1-c ~ 3e-4) that cancellation
//  costs ~3.5 of float32's ~7 digits, and the noise lands on everything
//  built from the curve. The rewritten form is a sum of non-negative
//  terms: fully conditioned. (1-c) itself is exact in float for c >= 0.5
//  (Sterbenz), i.e. precisely the sharp shapes that need it.

cheapCurveBase(c, x) = (log(c*x*x+(1-c))+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);

// --- 2c. curveScale: normalization factor ---
//
//  The difference between the curve at x=1 and x=0. Dividing by this
//  maps the raw curve to [0, 1] over the full phase range.

curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

// --- 2d. Normalized release and attack curves ---
//
//  cheapCurveRelease: starts at 0 (phase=0) and rises to 1 (phase=1).
//    This is the natural shape of the antiderivative, just normalized.
//
//  cheapCurveAttack: starts at 1 and falls to 0.
//    It's the horizontal+vertical flip of the release curve:
//      attack(x) = 1 - release(1 - x)
//
//  Why the flip? The same S-curve shape is reused for both directions.
//  For attack (envelope falling toward a transient), we want to start
//  fast and end slow ("ease out"). For release (envelope rising back up),
//  we want to start slow and end fast ("ease in"). Flipping gives us both
//  from a single curve definition.

cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x) = 1-cheapCurveRelease(c, 1-x);

// --- 2e. Derivatives of the curve ---
//
//  The derivative tells us the SPEED of the envelope at a given phase.
//  Since cheapCurveBase is the antiderivative of x(1-x)/(cx²+1-c),
//  the derivative of the normalized release curve is just that rational
//  function divided by curveScale.
//
//  derivativeBaseRelease:  unscaled (not divided by curveScale)
//  derivativeRelease:      scaled (the true derivative of the [0,1] curve)
//
//  The attack derivative uses the flip: d/dx attack(x) = release'(1-x),
//  but with the denominator adjusted for the substitution u = 1-x:
//    derivativeBaseAttack(c, x) = x(1-x) / (cx² + 1 - 2cx)

//  Denominators rewritten as sums of non-negative terms for float32
//  conditioning at sharp shapes (see the note above cheapCurveBase):
//    c*x*x+1-c       == c*x*x+(1-c)
//    c*x*x+1-2*c*x   == c*(x-1)*(x-1)+(1-c)
derivativeBaseRelease(c, x) = x*(1-x)/(c*x*x+(1-c));
derivativeBaseAttack(c, x) = x*(1-x)/(c*(x-1)*(x-1)+(1-c));

derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

// --- 2f. Peak phase of the derivative ---
//
//  The derivative is a bell-shaped hump (zero at phase 0 and 1, positive
//  in between). Its peak tells us where the envelope is moving fastest.
//  For attack, the peak is biased toward the start (fast initial plunge);
//  for release, it's biased toward the end.
//
//  These are used in the env loop to decide whether we're on the
//  "accelerating" or "decelerating" side of the curve — crucial for
//  picking the correct branch of the inverse derivative (see 2g).

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
peakPhaseRelease(c) = 1-peakPhaseAttack(c);
// mirror image

maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));

// --- 2g. Inverse derivatives ---
//
//  Given a target speed D, at what phase does the curve have that speed?
//
//  Because the derivative is a hump (rises then falls), there are TWO
//  solutions for any speed below the peak: one on the ascending ("bottom")
//  side and one on the descending ("top") side. The env loop must pick
//  the correct one based on whether the envelope is still accelerating
//  or already decelerating.
//
//  The math comes from solving  x(1-x)/(cx²+1-c) = D  for x, which is
//  a quadratic in x. The ± in the quadratic formula gives the two branches.
//
//  "Top" = larger phase (past the peak), "Bottom" = smaller phase (before peak).

inverseDerivativePart(c, D) = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));

inverseDerivativeTopRelease(c, D) = (1+inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeBottomRelease(c, D) = (1-inverseDerivativePart(c, D))/(2*(c*D+1));

// Attack inverses: horizontal flip of the release inverses
inverseDerivativeTopAttack(c, D) = 1-inverseDerivativeTopRelease(c, D);
inverseDerivativeBottomAttack(c, D) = 1-inverseDerivativeBottomRelease(c, D);

// ============================================================================
//  3. THE ENVELOPE FOLLOWER (env)
// ============================================================================
//
//  The heart of the algorithm. `env` is a feedback loop that moves its
//  output (`prev`) toward the target (`lookaheadX`) along the shaped curve.
//
//  Each sample, it:
//    a) Decides whether it's ATTACKING (prev > target, envelope falling)
//       or RELEASING (prev < target, envelope rising).
//    b) Computes the current PHASE (how far along the transition we are)
//       and the remaining SPAN (total distance of this transition).
//    c) Finds the phase whose curve-derivative matches the current speed
//       (velocity matching), so mid-transition target changes don't cause
//       discontinuities — the envelope seamlessly re-anchors onto the curve.
//    d) Steps phase forward by 1/N (where N = transition duration in samples)
//       and reads the curve's derivative at the new phase to get `speed`.
//    e) Computes delta = speed * step * totalStep, adds it to prev.
//
//  The velocity-matching trick (step c) is what makes this different from
//  a simple "play back the curve" approach. When the target moves mid-
//  transition, the envelope finds the point on the NEW curve that has the
//  same speed as it currently has, then continues from there. This gives
//  smooth, glitch-free transitions even with rapidly changing input.

// --- Smoother (raw; scaled att/rel/rel_hold are derived in Parameters) ---
attMs = SmootherGroup(hslider("[0]att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01));
attackShapeSlider = SmootherGroup(hslider("[1]attack shape", 0, 0, 1, 0.001));
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));

// Derived parameters
maxHold = 0.05;
// max release hold time (seconds)
maxSR = 48000;
maxHoldSamples = maxHold*maxSR;

att = attMs/1000;
rel = relMs/1000;
att_samples = att*ma.SR:max(1);
rel_samples = rel*ma.SR:max(1);

// Simplified process: test signal -> delay for latency comp, and the smoother
process = testSignal<:(_@int(att_samples), shapedSmoother(_));

shapedSmoother(x) = lookaheadX, delayedX:env~(_, _, _)
    with {
        // --- Lookahead via sliding minimum ---
        // Look att_samples into the future to find the minimum.
        // This lets the envelope start falling BEFORE a transient arrives.
        // Cost: att_samples of latency.
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHoldSamples);

        // The raw input delayed to align with the lookahead. Passed explicitly
        // to env so it doesn't become a free-variable (implicit extra input).
        delayedX = x@int(att_samples);

        // --- Precomputed per-shape constants (slider-rate, not audio-rate) ---
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
        //
        //  Inputs:
        //    lookaheadX    — the target value (from the sliding minimum)
        //    delayedX      — the raw input delayed by att_samples (hard ceiling)
        // =======================================================================
        env(prev, prevPhase, prevTotalStep, lookaheadX, delayedX) = result, newPhase, totalStep
            with {
                // --- Step 1: Are we attacking, releasing, or idle? ---
                //
                // attacking: envelope is above target, needs to fall
                // releasing: envelope is below target, needs to rise
                // If neither, the envelope is already at the target (idle).
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                // --- Step 2: Select shape parameters for current direction ---
                //
                // All these are slider-rate constants, precomputed above.
                // select2(releasing, attackVal, releaseVal) picks the right set.
                shape = select2(releasing, attackShape, releaseShape);
                invCurveScale = attackInvCurveScale+releasing*(releaseInvCurveScale-attackInvCurveScale);
                zeroVal = select2(releasing, attackZero, releaseZero);
                maxDerivVal = select2(releasing, attackMaxDerivBase, releaseMaxDerivBase);

                // Per-sample phase increment: 1/(duration_in_samples)
                step = select2(releasing, attStep, relStep);

                // --- Step 3: Compute totalStep (the full span of this transition) ---
                //
                // totalStep is the distance from the (fixed) start of the
                // current transition to the (moving) target. The envelope's
                // own progress lives in the phase, not here — so the only
                // thing that can change the span mid-transition is the target
                // moving. Integrate exactly that:
                //
                //   totalStep = prevTotalStep + (lookaheadX - lookaheadX')
                //
                // This is the old "hold" generalized: a static target gives a
                // zero increment and the span is held exactly; a moving
                // target shrinks or grows the span by exactly its own
                // movement, in either direction. Unlike reconstructing the
                // span through the curve (the old rawTotalStep =
                // remaining + prevTotalStep*fracDone), the state never
                // round-trips through cheapCurveBase, so there is no
                // accounting drift to re-inject — and no needsRecompute
                // heuristic to defeat: noise on the input just makes the
                // increments zero-mean jitter that sums to the target's net
                // movement by construction. (This also made fracDone and the
                // Step 3 curve evaluation dead code; they are gone.)
                //
                // On entry into a direction the span is simply the full
                // remaining distance. Entry is detected by the sign of
                // prevTotalStep: it carries the sign of the previous
                // direction (positive release, negative attack, zero idle).
                //
                // supportTotalStep is the smallest span whose curve can carry
                // the current speed (the curve derivative maxes out at
                // maxDerivVal). When the target jumps closer mid-transition
                // the span may shrink only down to this floor; below it,
                // clampedRatio would saturate and the envelope would slow
                // discontinuously. On the floor it instead decelerates along
                // the curve at the curve's steepest rate, and any resulting
                // overshoot is caught by the output clamps in Step 5. One
                // expression serves both directions: it floors the release
                // span (both positive) and ceils the attack span (both
                // negative).
                incrementalTotalStep = prevTotalStep+(lookaheadX-lookaheadX');
                entryTotalStep = lookaheadX-prev;
                supportTotalStep = prevSpeed/(maxDerivVal*invCurveScale*step);

                totalStep = select2(releasing,
                    min(select2(prevTotalStep>=0, incrementalTotalStep, entryTotalStep),
                        supportTotalStep),
                    max(select2(prevTotalStep<=0, incrementalTotalStep, entryTotalStep),
                        supportTotalStep))*active;

                // --- Step 4: Velocity matching ---
                //
                // What speed is the envelope currently moving at?
                // Use the ACTUAL output velocity, so that any external
                // algorithms that modify `prev` are automatically accounted for.
                prevSpeed = prev-prev';

                // Express current speed as a ratio of the new span
                // When idle, totalStep and prevSpeed are both 0.
                // Guard the denominator so speedRatio = 0/ε = 0
                // rather than 0/0 = NaN (which would persist in
                // prevPhase and poison all subsequent samples,
                // since IEEE 754 says NaN * 0 = NaN).
                speedRatio = prevSpeed/(totalStep*invCurveScale*step+(1-active)*1e-30);
                clampedRatio = max(0, min(speedRatio, maxDerivVal));

                // Find the phase on the curve that has this speed.
                // This is the KEY trick: instead of resetting the phase when the
                // target moves, we find where on the new curve our current speed
                // would naturally occur, then continue from there. Result: the
                // envelope seamlessly transitions to the new trajectory without
                // any discontinuity in velocity.
                //
                // Two branches (top/bottom) because the derivative curve is a hump:
                //   - "bottom" = ascending side (before peak, accelerating)
                //   - "top"    = descending side (after peak, decelerating)
                //
                // For ATTACK: gonnaMakeIt (will we overshoot?) picks the branch.
                //   If yes -> top (decel), if no -> bottom (accel).
                //
                // For RELEASE: gonnaMakeIt likewise picks the branch.
                gonnaDo(phase) = select2(releasing, cb, 1-cb)*totalStep
                    with {
                        cb = (cheapCurveBase(shape,
                            select2(releasing, 1-phase, phase))-zeroVal)*invCurveScale;
                    };

                // Euler-Maclaurin correction: the output is a right-endpoint
                // Riemann sum of the curve derivative (speed is read at
                // phaseAtMatchingSpeed+step), so the distance the discrete
                // loop will actually cover from a given phase falls short of
                // the continuous integral gonnaDo() computes by
                // ~(step/2)*derivNorm(phase)*totalStep — the boundary term of
                // the Euler-Maclaurin formula (the term at the far end
                // vanishes because the derivative is 0 at the curve endpoint).
                // Since totalStep carries the direction sign, one expression
                // covers both: it lowers the release projection and raises the
                // attack projection. derivativeBase at the inverse-derivative
                // phase IS clampedRatio by definition, so the correction is
                // free. Without it, `projected` drifts across lookaheadX along
                // the decelerating tail, gonnaMakeIt flips back to the
                // accelerating branch, and the envelope re-bursts: a kink on
                // release (caught by the delayedX clamp), overshoot below the
                // target on attack. Severe at sharp shapes, where one phase
                // step spans a large derivative range.
                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio)))+prev
                    -(0.5*step*clampedRatio*invCurveScale*totalStep);
                gonnaMakeIt = (projected>lookaheadX);

                phaseAtMatchingSpeed = select2(releasing,
                    // Attack branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        // accel
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    // decel
                    // Release branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        // accel
                        inverseDerivativeTopRelease(shape, clampedRatio)));
                // decel

                // --- Step 5: Advance phase and compute output ---
                //
                // newPhase = matched phase + one step, clamped to (step, 1-step)
                // to avoid the curve endpoints where the derivative is zero.
                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                // Read the curve's derivative at newPhase to get the current speed
                speed = totalStep*derivativeBaseRelease(shape,
                    select2(releasing, 1-newPhase, newPhase))*invCurveScale;

                // delta = speed * step = how much the envelope moves this sample
                delta = speed*step;

                // Final output: move toward target but never overshoot the raw
                // input (delayed by att_samples to align with the lookahead),
                // and never overshoot BELOW the attack target. Even with the
                // corrected projection, the attack arrives at near-peak speed
                // at sharp shapes (only ~a handful of samples of braking room
                // past the derivative peak), leaving a landing error of up to
                // one peak-speed sample; the lower clamp turns that into an
                // exact stop at lookaheadX. The bound min(prev, lookaheadX) is
                // inert outside of attack: while releasing or idle it equals
                // prev, which result never goes below. The two clamps cannot
                // conflict, since lookaheadX <= delayedX by construction
                // (the sliding-min window contains the delayed sample).
                result = min(prev+delta, delayedX):max(min(prev, lookaheadX));
            };
    };

// ============================================================================
//  GLOSSARY
// ============================================================================
//
//  phase:
//    A value in [0, 1] tracking progress through the current attack or
//    release transition. 0 = just started, 1 = done.
//
//  totalStep:
//    The full amplitude span (in linear gain) of the current transition,
//    from where it started to where it's heading. This scales the normalized
//    [0,1] curve to actual gain values. Maintained incrementally: held while
//    the target is static, moved by exactly (lookaheadX - lookaheadX') when
//    it isn't, recomputed as (lookaheadX - prev) on entry into a direction,
//    and never allowed to shrink past the span whose curve can still carry
//    the current speed (supportTotalStep).
//
//  step:
//    Phase increment per sample = 1 / (duration_in_samples). After
//    `duration` samples, phase has advanced from 0 to 1.
//
//  cheapCurveBase(c, x):
//    The raw antiderivative whose shape is controlled by c. Not normalized.
//
//  curveScale(c):
//    = cheapCurveBase(c,1) - cheapCurveBase(c,0). Normalizes the curve to [0,1].
//
//  derivativeBaseRelease(c, x):
//    The speed of the normalized release curve at phase x (unscaled by
//    curveScale — multiply by invCurveScale for the true derivative).
//
//  velocity matching:
//    When the target changes mid-transition, find the phase on the (possibly
//    new) curve whose derivative equals the envelope's current speed, then
//    continue from there. Ensures smooth, jerk-free transitions.
//
//  Euler-Maclaurin correction:
//    The per-sample deltas are a right-endpoint Riemann sum of the curve
//    derivative, which covers slightly less distance than the continuous
//    curve. `projected` subtracts the first-order boundary term
//    (step/2 * derivNorm * totalStep) so gonnaMakeIt compares against the
//    distance the discrete loop will actually cover; totalStep's sign makes
//    the same expression correct for both directions. Velocity matching pins
//    the integrator to right-endpoint Euler (the realized speed must
//    correspond exactly to the phase to resume from), so the correction must
//    live in the projection, not the integrator. The residual attack landing
//    error (switching granularity at near-peak speed) is caught by the
//    max(min(prev, lookaheadX)) clamp on the output.
//
//  gonnaMakeIt:
//    "If we take the accelerating branch of the inverse derivative, will
//    the resulting curve overshoot the target?" If yes, we should be on the
//    decelerating branch instead.
//
//  lookaheadX:
//    The sliding minimum of the input over the next att_samples. This is
//    the envelope's target: the lowest value it will need to reach within
//    the attack window.
//
//  AUC compensation:
//    (Omitted from this core file for clarity.) Adjusts the attack/release
//    duration so that changing the shape doesn't change the average level
//    of the envelope — sharper shapes get proportionally shorter durations.
