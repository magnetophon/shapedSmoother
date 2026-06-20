declare name "shapedSmoother_split";
declare version "0.16";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — attack and release shaped followers, fused core, v0.16
//
//  Derived from shapedSmoother_core v0.14. One shared machine, coreSmoother,
//  shapes a transition in either direction; the live direction selects the
//  curve parameters and the transcendentals run once. Two specializations
//  expose single-direction followers, and a both-directions instance is the
//  limiter:
//
//      attackSmoother(x, ceiling)  = coreSmoother(1, 0, ...)
//                           shapes DOWNWARD motion (attack curve),
//                           passes upward motion through instantly.
//      releaseSmoother(x, ceiling) = coreSmoother(0, 1, ...)
//                           shapes UPWARD motion (release curve),
//                           passes downward motion through instantly.
//      shapedSmoother              = coreSmoother(1, 1, ...) on slidingMin
//                           shapes BOTH directions (the composition below).
//
//  `ceiling` clamps the output to <= ceiling (the env <= dry limiter
//  guarantee, see RESTORED below); pass ma.INFINITY for no clamp.
//
//  The core is the v0.4 transition machine — trusted phase, velocity
//  matching, incremental/entry/support span, guarded landing, and the
//  support-overshoot clamp. A direction flip re-anchors by velocity match;
//  an inherited opposite-sign speed clamps that to a standstill (matchD = 0).
//
//  DROPPED relative to v0.14:
//
//    * THE BRAKE (v0.12-v0.14). Nothing decelerates a release into an
//      oncoming attack, so release-to-attack flips carry the inherited
//      speed and enter the attack from a standstill — the v0.4 velocity
//      corner (the v0.14 notes measured 600-1300 such corners per minute on
//      dense material). brakeEnable is gone: this file IS the brake-off end
//      of that dial. C1 across a flip is not maintained.
//
//  PRESERVED: the shaped curves, velocity matching, the reach rule, guarded
//  landing, float32 conditioning, and — through the slidingMin front-end in
//  the composition — anticipatory attack at the SAME latency (att_samples).
//  The dry delay realises the lookahead; the follower path adds none.
//
//  RESTORED: the limiter ceiling (envelope <= delayed dry input). v0.14's
//  in-loop min(_, delayedX) is back as a `ceiling` argument, fed back so a
//  biting clamp triggers velocity matching. The composition passes the
//  delayed dry as the ceiling. The per-target support-overshoot clamps are
//  kept regardless.
//
//  CPU: ONE curve evaluation per sample — 1 log + 1 atan + 1 sqrt, the same
//  as the unified v0.14 machine. The both-directions instance does NOT cost
//  two: the direction selects the scalar parameters and the transcendentals
//  run once (the shared-inversion idea from v0.7). Chaining two separate
//  followers (slidingMin -> releaseSmoother -> attackSmoother) instead would
//  evaluate two curves per sample — 2x. Verified against that chain: the
//  fused core matches it to float epsilon on steps, transients, ducks and
//  noise (where every handoff is at a standstill — the duck bottom, or a
//  release->attack flip), and to within ~0.05 only on a sustained gradual
//  ramp that reverses direction mid-glide, the one corner where a single
//  velocity-matched machine and two chained machines part ways. On limiter
//  material the two are interchangeable; the fused core is that chain,
//  collapsed.
//
//  The composition at the bottom is
//      coreSmoother(1, 1) on slidingMin(att), ceiling = delayed dry
//  giving anticipatory attack-shaped ducks and release-shaped recovery,
//  with the flip corner the brake used to remove.
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
//  coreSmoother — shaped follower, one machine for both directions
//
//    shapeDown : 1 = shape DOWNWARD motion (attack curve); 0 = instant down
//    shapeUp   : 1 = shape UPWARD   motion (release curve); 0 = instant up
//
//  The live direction (prev > target moves down/attack, prev < target moves
//  up/release) selects the SCALAR curve parameters, and the transcendentals
//  are applied once to the selected scalars. So a both-directions instance
//  (shapeDown = shapeUp = 1) costs ONE log + atan + sqrt per sample — the
//  same as the unified machine, not two. A single-direction instance has
//  isRel constant (folded):
//
//    attackSmoother  = coreSmoother(1, 0, ...)   shape down, instant up
//    releaseSmoother = coreSmoother(0, 1, ...)   shape up,   instant down
//    the limiter     = coreSmoother(1, 1, ...)   shape both  (shapedSmoother)
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
coreSmoother(shapeDown, shapeUp, target, ceiling) = (target, ceiling) : loop ~ si.bus(4) : (_, !, !, !)
with {
    loop(prev, prevPhase, prevTotalStep, prevExpected, target, ceiling)
        = result, newPhase, totalStepOut, expected
    with {
        prevSpeed = prev-prev';

        // Per-sample direction. attacking/releasing gate on the shape mask
        // and the rate, so a masked-off direction never "moves" — it passes
        // through instantly (result clamp below).
        attacking = (prev>target) & (att>0) & (shapeDown!=0);
        releasing = (prev<target) & (rel>0) & (shapeUp!=0);
        moving = attacking | releasing;
        active = moving;

        // Curve-parameter selector. Single-direction: the constant shapeUp
        // (0 attack, 1 release), folded. Both-directions: the live direction.
        dualMode = (shapeDown!=0) & (shapeUp!=0);
        isRel = select2(dualMode, shapeUp, prev<target);

        // Direction-selected scalars (select on slider-rate constants; the
        // transcendentals below then run once on the chosen values).
        shp       = select2(isRel, attackShape, releaseShape);
        invCS     = select2(isRel, attackInvCurveScale, releaseInvCurveScale);
        zeroV     = select2(isRel, attackZero, releaseZero);
        maxDB     = select2(isRel, attackMaxDerivBase, releaseMaxDerivBase);
        stp       = select2(isRel, attStep, relStep);
        suppRecip = select2(isRel, attackSupportRecip, releaseSupportRecip);
        kV        = select2(isRel, sqrt((1/attackShape)-1), sqrt((1/releaseShape)-1));
        twoC      = select2(isRel, 2*attackShape, 2*releaseShape);
        oneMinusC = select2(isRel, 1-attackShape, 1-releaseShape);
        halfStep  = 0.5*stp;

        // --- Curve helpers (ONE direction's transcendentals per sample) ---
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
        sameDirection = (attacking==attacking') & (releasing==releasing');
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
        // The final :min(ceiling) restores the v0.14 limiter guarantee
        // (envelope <= delayed dry input). It is fed back as prev, so a
        // biting ceiling disturbs the output and triggers velocity
        // matching next sample, exactly like v0.14's in-loop min(_,
        // delayedX). ceiling = ma.INFINITY makes it a no-op (the bare
        // attackSmoother/releaseSmoother); the composition feeds the
        // delayed dry to the final attack stage.
        preClamp = select2(landed, expected, target);
        movingVal = select2(isRel,
            (preClamp:min(prev):max(target)),
            (preClamp:max(prev):min(target)));
        result = select2(moving, target, movingVal) : min(ceiling);

        // State resets on landing OR while passing through instantly
        // (active = 0 zeroes totalStep and newPhaseRaw).
        newPhase = newPhaseRaw*(1-phaseDone);
        totalStepOut = totalStep*(1-phaseDone);
    };
};

// The two public followers. `ceiling` clamps the output to <= ceiling
// (the env <= dry limiter guarantee), fed back so a biting clamp triggers
// velocity matching next sample. For an unclamped follower pass
// ceiling = ma.INFINITY.
// The two single-direction followers, as specializations of coreSmoother.
// `ceiling` clamps the output to <= ceiling (the env <= dry limiter
// guarantee), fed back so a biting clamp triggers velocity matching next
// sample; pass ma.INFINITY for no clamp. With one shape flag off, isRel is
// constant and the unused direction's curve folds away (1 transcendental).
attackSmoother(x, ceiling)  = coreSmoother(1, 0, x, ceiling);
releaseSmoother(x, ceiling) = coreSmoother(0, 1, x, ceiling);

// ============================================================================
//  COMPOSITION
//
//  slidingMin gives the anticipatory target (it drops the moment a new
//  minimum enters the window, att_samples before the transient reaches the
//  delayed dry). A single both-directions coreSmoother then shapes the duck
//  DOWN onto it (attack curve) and the recovery UP off it (release curve).
//  Latency = lookaheadSamples, spent on the dry delay, not the follower.
//
//  Why fused, not slidingMin -> releaseSmoother -> attackSmoother:
//  two separate followers each evaluate their own curve every sample, and
//  Faust computes both even though, here, only one direction is ever moving
//  at a time (the duck completes before the recovery starts; a transient
//  during a recovery flips release->attack from a standstill). That is a
//  2x transcendental bill for almost no behavioural difference. The fused
//  follower picks the live direction's parameters and evaluates the curve
//  ONCE per sample, so it costs 1 log + atan + sqrt instead of 2 — and is
//  equivalent to the release->attack chain to float epsilon wherever the
//  handoff is at a standstill (the duck bottom at ~zero velocity; a
//  release->attack flip, the v0.4 corner), i.e. on all transient/step/noise
//  material. Only a sustained gradual ramp reversing mid-glide differs, and
//  only slightly (~0.05). attackSmoother and
//  releaseSmoother remain above for genuinely independent single-direction
//  use; chaining them is the slow path.
//
//  env <= dry: the fused follower is given the delayed dry as its ceiling,
//  so the envelope never exceeds the aligned input. It is fed back, so a
//  biting clamp trips velocity matching next sample (the v0.14 guarantee).
// ============================================================================
shapedSmoother(x) = coreSmoother(1, 1, mn, dry)   // shape both directions, ceiling = delayed dry
with {
    dry = x @ lookaheadSamples;
    mn  = x : ba.slidingMin(int(att_samples)+1, 1+maxLookaheadSamples);
};

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
