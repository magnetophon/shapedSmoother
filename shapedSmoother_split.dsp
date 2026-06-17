declare name "shapedSmoother_split";
declare version "0.15";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — attack and release as two independent followers, v0.15
//
//  Derived from shapedSmoother_core v0.14 by splitting the unified
//  envelope state machine into two single-direction shaped followers
//  that compose:
//
//      attackSmoother(x)  — shapes DOWNWARD motion (attack curve),
//                           passes upward motion through instantly.
//      releaseSmoother(x) — shapes UPWARD motion (release curve),
//                           passes downward motion through instantly.
//
//  Each is the v0.4 transition machine — trusted phase, velocity
//  matching, incremental/entry/support span, guarded landing, and the
//  support-overshoot clamp — with the direction PINNED and the opposite
//  direction passed through in a single sample. There is NO shared state
//  between the two: they are genuinely independent and can be reordered,
//  reused, or replaced separately.
//
//  DROPPED relative to v0.14 (the price of independence):
//
//    * THE BRAKE (v0.12-v0.14). A release no longer sees the upcoming
//      attack, so it cannot decelerate into it. Release-to-attack flips
//      therefore carry the inherited speed — the v0.4 velocity corner is
//      back (the v0.14 notes measured 600-1300 such corners per minute on
//      dense material, up to 77% of peak speed). brakeEnable is gone:
//      this file IS the brake-off end of that dial, structurally.
//
//    * CROSS-DIRECTION VELOCITY MATCHING. Each follower velocity-matches
//      only within its own direction. At a flip the downstream follower
//      simply passes the upstream motion through; the new direction
//      enters from a standstill, because an inherited opposite-sign speed
//      clamps the velocity match to phase 0 (matchD = 0). C1 across a
//      flip is not maintained — that was the brake's job.
//
//    * THE LIMITER CEILING (output <= delayed dry input). v0.14's
//      min(_, delayedX) coupled the envelope to the dry signal. The
//      followers clamp only to their OWN target — the support-overshoot
//      catch, which IS required and is kept. To restore env <= dry, give
//      attackSmoother the delayed dry as an extra ceiling input and
//      append :min(ceiling) to the attack branch of its result clamp;
//      see the note there.
//
//  PRESERVED: the shaped curves, per-direction velocity matching, the
//  reach rule, guarded landing, float32 conditioning, and — through the
//  slidingMin front-end in the shapedSmoother composition below —
//  anticipatory attack at the SAME latency (att_samples). The dry delay
//  realises the lookahead; the follower path adds none.
//
//  CPU: two curve evaluations per sample instead of v0.14's shared one,
//  i.e. 2 log + 2 atan + 2 sqrt where the unified machine had 1/1/1. That
//  doubling is the unavoidable cost of two independent machines and is
//  exactly what the shared-inversion work (v0.7) existed to avoid.
//
//  The composition at the bottom wires
//      slidingMin(att) -> attackSmoother -> releaseSmoother
//  reproducing the v0.4-family behaviour: anticipatory attack-shaped
//  ducks, release-shaped recovery, with the flip corner the brake used to
//  remove.
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

testSignal1 = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };
testSignal2 = os.lf_squarewave(testFreq)*0.5;

// ============================================================================
//  THE SHAPED CURVES  (unchanged from v0.14)
// ============================================================================

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));

cheapCurveBase(c, x) = (log(c*x*x+(1-c))+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);

curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

derivativeBaseRelease(c, x) = x*(1-x)/(c*x*x+(1-c));
derivativeBaseAttack(c, x) = x*(1-x)/(c*(x-1)*(x-1)+(1-c));

derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);

maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));

inverseDerivativePart(c, D) = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));

// ============================================================================
//  CONTROLS AND DERIVED PARAMETERS
// ============================================================================

attMs = SmootherGroup(hslider("[0]att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01));
attackShapeSlider = SmootherGroup(hslider("[1]attack shape", 0, 0, 1, 0.001));
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));

maxHold = 0.05;
maxSR = 48000;
maxHoldSamples = maxHold*maxSR;
maxLookaheadSamples = 2*maxHoldSamples+1;

att = attMs/1000;
rel = relMs/1000;
att_samples = att*ma.SR:max(1);
rel_samples = rel*ma.SR:max(1);

// Latency = the attack window. No brake window: an independent release
// has nothing to brake against (see header). The follower signal path
// adds no delay; this delay aligns the dry signal so the anticipatory
// duck lands exactly when its transient arrives.
lookaheadSamples = min(int(att_samples), int(maxLookaheadSamples));

// --- Per-shape constants (slider-rate), shared by both followers ---
attackShape = shapeMap(attackShapeSlider);
releaseShape = shapeMap(releaseShapeSlider);

attackInvCurveScale = 1/curveScale(attackShape);
releaseInvCurveScale = 1/curveScale(releaseShape);

attackZero = cheapCurveBase(attackShape, 0);
releaseZero = cheapCurveBase(releaseShape, 0);

attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

// Fencepost: a transition occupies the CLOSED interval of N+1 output
// samples, so the phase advances in N+1 steps of 1/(N+1). (Same reasoning
// as v0.4; it is what makes the duck land exactly when the delayed dry
// transient arrives.)
attStep = 1/((att*ma.SR):max(1)+1);
relStep = 1/((rel*ma.SR):max(1)+1);

// Reciprocal of the support-floor denominator, hoisted to slider rate so
// the loop multiplies instead of dividing on the recurrence path (v0.7).
attackSupportRecip = 1/(attackMaxDerivBase*attackInvCurveScale*attStep);
releaseSupportRecip = 1/(releaseMaxDerivBase*releaseInvCurveScale*relStep);

// ============================================================================
//  dirSmoother — one-directional shaped follower
//
//    isRel = 0 : ATTACK  — shapes downward motion, instant up.
//    isRel = 1 : RELEASE — shapes upward motion, instant down.
//
//  With isRel a literal, every select2(isRel, ...) below constant-folds,
//  so each instance evaluates exactly one direction's curve.
//
//  State carried between samples:
//    prev          — current envelope value.
//    prevPhase     — position on the curve, [0, 1].
//    prevTotalStep — span of the current transition; == 0 exactly iff
//                    idle (the completion reset writes literal 0).
//    prevExpected  — the value this loop computed last sample; if
//                    prev != prevExpected the output was clamped or
//                    modified outside, so velocity-match instead of
//                    trusting the stored phase.
// ============================================================================
dirSmoother(isRel, target) = target : loop ~ si.bus(4) : (_, !, !, !)
with {
    // Direction-fixed per-shape selections (constant-folded per instance).
    shp       = select2(isRel, attackShape, releaseShape);
    invCS     = select2(isRel, attackInvCurveScale, releaseInvCurveScale);
    zeroV     = select2(isRel, attackZero, releaseZero);
    maxDB     = select2(isRel, attackMaxDerivBase, releaseMaxDerivBase);
    stp       = select2(isRel, attStep, relStep);
    suppRecip = select2(isRel, attackSupportRecip, releaseSupportRecip);
    kV        = select2(isRel, sqrt((1/attackShape)-1), sqrt((1/releaseShape)-1));
    twoC      = select2(isRel, 2*attackShape, 2*releaseShape);
    oneMinusC = select2(isRel, 1-attackShape, 1-releaseShape);
    rateOK    = select2(isRel, att>0, rel>0);
    halfStep  = 0.5*stp;

    loop(prev, prevPhase, prevTotalStep, prevExpected, target)
        = result, newPhase, totalStepOut, expected
    with {
        prevSpeed = prev-prev';

        // Smoothed direction active?  attack: prev>target ; release:
        // prev<target. The opposite direction is the instant branch
        // (result, below).
        moving = select2(isRel, prev>target, prev<target) & rateOK;
        active = moving;

        // --- Curve helpers (one direction's transcendentals) ---
        cbAt(p) = ((log(shp*xx*xx+oneMinusC)
            +2*kV*atan(xx/kV)-2*xx)/twoC-zeroV)*invCS
        with { xx = select2(isRel, 1-p, p); };
        fdOf(p) = select2(isRel, 1-cbAt(p), cbAt(p));
        derivOf(p) = derivativeBaseRelease(shp, select2(isRel, 1-p, p))*invCS;

        // --- Step 1: totalStep (incremental; entry on restart; support) ---
        // No brake: the target is the plain input. The support floor uses
        // the SIGNED prevSpeed on BOTH sides (v0.4): an inherited
        // opposite-sign speed must not inflate the span — velocity
        // matching clamps such an entry to a standstill anyway, and the
        // v0.12 two-sided floor only existed for the (now absent)
        // turnaround.
        incrementalTotalStep = prevTotalStep+(target-target');
        entryTotalStep = target-prev;
        supportTotalStep = prevSpeed*suppRecip;
        totalStep = select2(isRel,
            min(select2(prevTotalStep>=0, incrementalTotalStep, entryTotalStep),
                supportTotalStep),
            max(select2(prevTotalStep<=0, incrementalTotalStep, entryTotalStep),
                supportTotalStep))*active;

        // --- Step 2: trusted phase vs velocity matching ---
        sameTarget = (target==target');
        sameDirection = (moving==moving');
        undisturbed = (prev==prevExpected);
        trusted = sameTarget&sameDirection&undisturbed&(prevTotalStep!=0);

        // --- Step 3: velocity matching (shared inversion, v0.7) ---
        matchCorr = halfStep*(sameDirection&(prevTotalStep!=0));
        speedRatio = prevSpeed/(totalStep*invCS*stp+(1-active)*1e-30);

        matchD = max(0, min(speedRatio, maxDB));
        partD = inverseDerivativePart(shp, matchD);
        denD = 2*(shp*matchD+1);
        topReleaseD = (1+partD)/denD;
        bottomReleaseD = (1-partD)/denD;

        lateAnchor = select2(isRel, 1-bottomReleaseD, topReleaseD)+matchCorr:min(1);
        gonnaDo(phase) = (1-fdOf(phase))*totalStep;
        projected = gonnaDo(lateAnchor)+prev;
        gonnaMakeIt = projected>target;

        // Reach rule (v0.4): pick the accelerating anchor unless the
        // decelerating anchor still reaches the target (overshoot is
        // clampable; stranding short arrives late).
        forwardPos = select2(isRel,
            select2(gonnaMakeIt, 1-bottomReleaseD, 1-topReleaseD),
            select2(gonnaMakeIt, bottomReleaseD, topReleaseD))
            +matchCorr:max(0):min(1-stp);
        matchedPos = forwardPos;

        // --- Step 4: advance phase ---
        anchor = select2(trusted, matchedPos, (prevPhase:min(1-stp)));
        newPhaseRaw = (anchor+stp)*active;

        // --- Step 5: curve increment (2-pt Gauss-Legendre); no brake fade ---
        glA = 0.21132486540518713;
        glB = 0.78867513459481287;
        delta = totalStep*stp*0.5
            *(derivOf(anchor+glA*stp)+derivOf(anchor+glB*stp));
        expected = prev+delta;

        // --- Step 6: guarded landing ---
        phaseDone = active&(newPhaseRaw>=(1-halfStep));
        guardSize = abs(totalStep)*maxDB*invCS*stp;
        landed = phaseDone&(abs(target-expected)<=guardSize);

        // --- Step 7: output ---
        // Smoothed direction: clamp to [target, prev] (attack) or
        // [prev, target] (release). The target side is the
        // support-overshoot catcher; when it fires the output differs from
        // expected, so the next sample velocity-matches. Opposite
        // direction: jump to target in a single sample (instant).
        //
        // To restore env <= dry (the dropped limiter ceiling), give this
        // function the delayed dry as a ceiling input and append
        // :min(ceiling) to the attack branch below.
        preClamp = select2(landed, expected, target);
        movingVal = select2(isRel,
            (preClamp:min(prev):max(target)),
            (preClamp:max(prev):min(target)));
        result = select2(moving, target, movingVal);

        // State resets on landing OR while passing through instantly
        // (active = 0 zeroes totalStep and newPhaseRaw).
        newPhase = newPhaseRaw*(1-phaseDone);
        totalStepOut = totalStep*(1-phaseDone);
    };
};

attackSmoother(x)  = dirSmoother(0, x);
releaseSmoother(x) = dirSmoother(1, x);

// ============================================================================
//  COMPOSITION
//
//  slidingMin gives the anticipatory attack target (it drops the moment a
//  new minimum enters the window, att_samples before the transient reaches
//  the delayed dry). attackSmoother shapes the duck onto it; releaseSmoother
//  shapes the recovery. Latency = lookaheadSamples, spent on the dry delay,
//  not the follower path.
//
//  Swap the order, drop a stage, or feed either follower a different
//  target — they are independent.
// ============================================================================
shapedSmoother(x) = x
    : ba.slidingMin(int(att_samples)+1, 1+maxLookaheadSamples)
    : attackSmoother
    : releaseSmoother;

process = testSignal<:(_@lookaheadSamples, shapedSmoother(_));

// ============================================================================
//  GLOSSARY  (terms specific to the brake/turnaround are gone with them)
// ============================================================================
//
//  phase:        Progress through the current transition, [0, 1]; advances
//                by step per sample, trusted between discontinuities and
//                velocity-matched across them.
//  totalStep:    Full signed span of the current leg; held while the target
//                is static, moved by (target - target') when it isn't,
//                (re)entered as (target - prev), floored by the span whose
//                curve can still carry the current speed (supportTotalStep).
//  step:         Phase increment = 1 / (duration_in_samples).
//  fdOf(p):      Fraction of the leg completed at phase p; position =
//                start + totalStep*fdOf(p). The delta integrates fdOf' over
//                one step with 2-pt Gauss-Legendre.
//  velocity matching: Across a discontinuity, find the phase whose
//                derivative equals the measured speed and continue from
//                there. Opposite-sign speed clamps to phase 0 (the leg
//                start) — e.g. an attack entered from an upward pass-through.
//  gonnaMakeIt:  "Would the decelerating branch still reach the target?"
//                Picks decelerating-if-it-reaches, else accelerating.
//  landing:      On phase completion the leg resets; a forward leg snaps to
//                the target only within one peak-speed sample, else the
//                leftover gets a fresh velocity-matched leg.
// ============================================================================
