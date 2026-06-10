declare name "shapedSmoother_core";
declare version "0.5";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — core algorithm, v0.5
//
//  CHANGES vs v0.4: THE TURNAROUND.
//
//  An attack arriving while a release was still in flight hit the
//  clampedRatio floor: the opposite-sign speed matched to phase 0, the
//  envelope stopped dead in one sample and re-accelerated from rest — an
//  instant change of direction. v0.5 buys extra lookahead and spends it on
//  turning around smoothly instead.
//
//  1. MIRRORED ATTACK CURVE / NEGATIVE PHASE. The attack's fraction-done
//     function is extended to negative phase by even reflection,
//     fd(p) := fd(|p|), so position(p) = apex + totalStep*fd(|p|), and the
//     chain rule flips the velocity's sign across zero: the derivative of
//     the extension is signum(p)*fd'(|p|). Advancing the phase from -q
//     through 0 to 1 then traces ONE trajectory: rise, decelerate to a
//     standstill at the apex (phase 0), then the unmodified attack down to
//     the target. The velocity is continuous everywhere, including the
//     apex (zero from both sides; fd is locally parabolic there, so the
//     2-point quadrature crosses it without noticing). Cost: abs() on the
//     curve evaluations, signum() on the derivative.
//
//  2. ARRIVAL COUNTDOWN. A sliding minimum can only drop because the
//     ENTERING sample is the window's new minimum, and the entering sample
//     reaches delayedX in exactly lookaheadSamples. So one counter, armed
//     when a drop lands below the envelope while nothing is pending,
//     counts down to the moment the pending minimum must be fully ducked.
//     Further drops while a countdown is pending do NOT re-arm it: the
//     FIRST deadline is the binding one, and mid-transition retargets
//     already land early (safe), exactly as in v0.4.
//
//  3. JUST-IN-TIME FLIP. Wanting to attack no longer means attacking NOW.
//     The flip is held while
//         arrival > att_samples + matchedDepth*(att_samples+1),
//     i.e. while there is still time to start the turnaround later. During
//     the hold the release simply continues — its span frozen against the
//     (already dropped) target, its phase trusted — and if it completes it
//     parks: the landing guard refuses the snap and the state resets to
//     idle. On the flip sample the mirror anchor is placed at the GRID
//     depth, -(arrival-att_samples)*attStep, which makes the landing
//     coincide with the arrival BY CONSTRUCTION: the traversal from -q to
//     1 takes att_samples + q/attStep samples, which is exactly what is
//     left on the countdown when the gate trips. Slow releases flip late
//     with a shallow mirror; an idle attack degenerates to depth 0 and
//     flips with exactly att_samples left — v0.4's timing. At
//     turnaround = 0 the whole feature is inert and v0.4's behavior
//     (corner included) returns: the A/B is one slider.
//
//  4. ENTRY SPAN CORRECTION. The mirror leg retraces ground: anchored at
//     -q, the total forward travel telescopes to totalStep*(1-fd(q)), not
//     totalStep. A turnaround entry therefore divides the entry span by
//     (1-fd(q)), so the trajectory lands ON the target instead of fd(q)
//     short of it. q is the grid depth, known before the span — no
//     circularity. The gate's depth estimate still uses the UNcorrected
//     span, so extreme turnarounds (incoming speed near the attack curve's
//     peak, soft shapes) carry a bounded slope error at the flip instead
//     of v0.4's full stop; if field testing shows it, one fixed-point
//     refinement of the depth (re-match against the corrected span) is the
//     known fix, at the price of a second curve evaluation.
//
//  5. SENTINEL MIGRATION. Phase 0 is now a legitimate mid-transition value
//     (the apex), so "no transition in progress" moves from
//     prevPhase == 0 to prevTotalStep == 0 — a value the completion reset
//     already establishes and which no live transition holds. The trusted
//     gate and the matchCorr gate use the new sentinel; the apex sample
//     stays trusted and the phase walks through 0 unmolested. (The phase
//     hits 0.0 exactly there: it is accumulated in whole steps from a
//     whole-step anchor, and -step+step is exact in float.)
//
//  6. TWO-SIDED SUPPORT FLOOR (attack). During the mirror leg the speed
//     has the opposite sign to the span, so the attack-side floor uses
//     -abs(prevSpeed): a span born in a turnaround can always carry the
//     speed it inherits. The release side keeps the signed prevSpeed —
//     entry from an attack must stay inert there.
//
//  COSTS. Latency rises from att_samples to
//  att_samples + turnaround*(att_samples+1): doubled lookahead at the
//  default turnaround = 1. Per sample it adds one inverse-derivative (the
//  gate), one curve evaluation (the span correction), a counter, and
//  abs/signum — roughly the budget v0.4's note B freed up. The sliding-min
//  tree allocates for the doubled window.
//
//  Everything else — cheap 2-pt Gauss-Legendre increments, trusted phase,
//  incremental span, guarded landing, float32 conditioning, the fencepost
//  and landing-slack reasoning — carries over from v0.4 unchanged.
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
turnaroundSlider = SmootherGroup(hslider("[2]turnaround", 1, 0, 1, 0.001));
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));

// Derived parameters
maxHold = 0.05;
maxSR = 48000;
maxHoldSamples = maxHold*maxSR;
maxLookaheadSamples = 2*maxHoldSamples+1;

att = attMs/1000;
rel = relMs/1000;
att_samples = att*ma.SR:max(1);
rel_samples = rel*ma.SR:max(1);

// Turnaround budget. The deepest mirror anchor velocity matching can ever
// produce is peakPhaseAttack < 1, whose traversal costs strictly less than
// one extra attack duration; extraSamples = att_samples+1 (turnaround = 1)
// therefore covers every matched turnaround for every shape, without tying
// the latency to the shape slider. lookaheadSamples is the single source
// of truth for the window length, the delay, AND the countdown reset: the
// drop<->arrival bookkeeping in Step 0 is only exact if all three agree.
extraSamples = int(turnaroundSlider*(att_samples+1));
lookaheadSamples = int(att_samples)+extraSamples;

process = testSignal<:(_@lookaheadSamples, shapedSmoother(_));

shapedSmoother(x) = lookaheadX, delayedX:env~(_, _, _, _, _):(_, _, _, !, !)
    with {
        lookaheadX = x:ba.slidingMin(lookaheadSamples+1, 1+maxLookaheadSamples);
        delayedX = x@lookaheadSamples;

        // --- Precomputed per-shape constants (slider-rate) ---
        attackShape = shapeMap(attackShapeSlider);
        releaseShape = shapeMap(releaseShapeSlider);

        attackInvCurveScale = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackZero = cheapCurveBase(attackShape, 0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
        releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

        attackPeakPhase = peakPhaseAttack(attackShape);

        // Fencepost: a transition occupies the CLOSED interval of N+1 output
        // samples [entry .. entry+N], so the phase advances in N+1 steps of
        // 1/(N+1). With 1/N the attack landed one sample early: lookaheadX
        // drops at m (the new minimum enters the window immediately) while
        // delayedX delivers the transient at m+480, and N increments
        // starting ON the entry sample complete at m+479. With 1/(N+1) the
        // attack's landing snap fires on exactly the sample where delayedX
        // presents the transient, and a release elapses exactly
        // rel_samples from onset to target. (In v0.5 the same arithmetic is
        // what makes the grid anchor of Step 5 land exactly when the
        // arrival countdown of Step 0 reaches zero.)
        attStep = 1/((att*ma.SR):max(1)+1);
        relStep = 1/((rel*ma.SR):max(1)+1);

        // =======================================================================
        //  env: the recursive core
        //
        //  State carried between samples:
        //    prev          — current envelope value
        //    prevPhase     — where we are on the curve. [0, 1] for a release;
        //                    [-peakPhaseAttack, 1] for an attack, where the
        //                    negative leg is the mirrored deceleration of a
        //                    turnaround and 0 is its apex.
        //    prevTotalStep — total span of the current transition. Also the
        //                    idle sentinel: == 0 exactly iff no transition is
        //                    in progress (the completion reset establishes it,
        //                    and phase 0 can no longer serve — it is the apex).
        //    prevExpected  — the value THIS loop computed last sample.
        //                    prev != prevExpected means something outside the
        //                    loop (the output clamps, or an external
        //                    algorithm) modified the envelope, so the stored
        //                    phase can no longer be trusted and we must
        //                    velocity-match.
        //    prevArrival   — samples until the pending window minimum reaches
        //                    delayedX; <= 0 when nothing is pending.
        // =======================================================================
        env(prev, prevPhase, prevTotalStep, prevExpected, prevArrival, lookaheadX, delayedX) = result, newPhase, totalStepOut, expected, arrival
            with {
                prevSpeed = prev-prev';

                // --- Step 0: arrival countdown ---
                //
                // A sliding minimum drops on a sample only because the
                // ENTERING input sample is the window's new minimum (the
                // window gains one sample and loses one; the min cannot
                // drop by losing a sample). The entering sample reaches
                // delayedX after exactly lookaheadSamples. So: on a drop
                // that lands below the envelope, arm the countdown to
                // lookaheadSamples; it then names the sample on which the
                // envelope must have fully ducked.
                //
                // Drops while a countdown is pending do not re-arm it. The
                // first deadline is the binding one; later, deeper minima
                // ride along as mid-transition retargets and land early
                // (safe), exactly as in v0.4. Drops that stay above the
                // envelope bind nothing (the envelope is already below
                // them) and do not arm.
                drop = lookaheadX<lookaheadX';
                arm = drop&(lookaheadX<prev)&(prevArrival<=0);
                arrival = select2(arm, max(prevArrival-1, 0), lookaheadSamples);

                // --- Step 1: direction, gated by the turnaround schedule ---
                //
                // Wanting to attack no longer means attacking NOW. needDepth
                // is the mirror anchor that velocity matching would pick for
                // the CURRENT speed (zero for an idle entry, or any speed
                // already pointing down); its traversal plus the attack
                // proper costs
                //     needSamples = att_samples + needDepth*(att_samples+1).
                // While the countdown still exceeds that, the attack is
                // HELD: flipping later both lets the release run longer and
                // is what makes the landing exact. The gate trips at
                // arrival <= needSamples; with needDepth = 0 that is
                // arrival == att_samples — v0.4's start time to the sample.
                //
                // During the hold the in-flight release (if any) continues:
                // span frozen (Step 3), phase trusted (Step 4). If it
                // completes first, the landing guard refuses to snap across
                // the dropped target and the state parks at idle until the
                // gate trips.
                wantsAttack = (prev>lookaheadX)&(att>0);
                needRatio = prevSpeed/(max(prev-lookaheadX, 1e-30)*attackInvCurveScale*attStep);
                needDepth = inverseDerivativeTopAttack(attackShape,
                    max(0, min(needRatio, attackMaxDerivBase)));
                needSamples = att_samples+needDepth*(att_samples+1);
                held = wantsAttack&(arrival>needSamples);

                attacking = wantsAttack&(1-held);
                releasing = ((prev<lookaheadX)|(held&(prevTotalStep>0)))&(rel>0);
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
                //          For negative (mirror) phases the position uses
                //          fdOf(|p|); only the derivative needs the sign
                //          (see derivM in Step 7).
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
                // movement. During a HOLD the increment is gated off: the
                // held release's target has already dropped away, and it
                // should finish (or be flipped out of) the trajectory it was
                // on, not chase a target that now belongs to the attack.
                //
                // On entry into a direction the span is the full remaining
                // distance — divided, for a turnaround entry, by
                // (1 - fd(grid depth)): the mirror leg retraces ground, so
                // the forward travel from anchor -q telescopes to
                // totalStep*(1-fd(q)). Without the correction the landing
                // would fall fd(q)*totalStep short and the guard would have
                // to re-transition the leftover — a late duck. q is the grid
                // depth (known before the span; no circularity), capped at
                // peakPhaseAttack both because matching never anchors deeper
                // and because the cap bounds the divisor away from 0.
                //
                // supportTotalStep is the smallest span whose curve can
                // carry the current speed. The attack side uses
                // -abs(prevSpeed): during the mirror leg the speed has the
                // opposite sign to the span, and a turnaround must be born
                // with a span that can carry the speed it inherits. The
                // release side keeps the signed prevSpeed so that entry
                // from an attack stays inert, as in v0.4. On the floor the
                // envelope decelerates along the curve at the curve's
                // steepest rate; the deliberate overshoot this allows is
                // caught by the output clamps in Step 7.
                mirrorBudget = max(arrival-att_samples, 0)*attStep;
                corrBudget = min(mirrorBudget, attackPeakPhase);
                turnEntry = wantsAttack&(prevSpeed>0);
                entryScale = 1-cheapCurveAttack(attackShape, corrBudget)*turnEntry;

                incrementalTotalStep = prevTotalStep+(lookaheadX-lookaheadX')*(1-held);
                entryTotalStep = (lookaheadX-prev)/entryScale;
                supportTotalStep = select2(releasing, 0-abs(prevSpeed), prevSpeed)
                    /(maxDerivVal*invCurveScale*step);

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
                //     unless we are HOLDING, where target motion belongs to
                //     the upcoming attack and the held release ignores it,
                //   - the direction flipped,             (sameDirection fails)
                //   - the envelope was modified outside  (undisturbed fails)
                //     this loop (external algorithms, or our own output
                //     clamps having altered the previous output),
                //   - or no transition was in progress.  (prevTotalStep == 0)
                // The last test moved off prevPhase == 0: the apex of a
                // turnaround sits at phase 0.0 exactly, and it is a trusted
                // mid-transition sample, not an entry.
                sameTarget = (lookaheadX==lookaheadX');
                sameDirection = (releasing==releasing')&(attacking==attacking');
                undisturbed = (prev==prevExpected);
                trusted = (sameTarget|held)&sameDirection&undisturbed&(prevTotalStep!=0);

                // --- Step 5: velocity matching (for untrusted samples) ---
                //
                // Since the emitted delta is an exact curve increment, the
                // measured speed is the MEAN derivative over the last step.
                // The inverse-derivative functions invert the POINT
                // derivative, so they return (approximately) the midpoint of
                // that step. Adding halfStep converts the recovered midpoint
                // into the current position — on either side of the apex:
                // the mirror phase advances by +step like any other, so the
                // recovered |midpoint| is current depth + halfStep and the
                // correction is additive after negation too.
                //
                // The correction only makes sense when a previous interval
                // in the same direction actually exists; the gate is the
                // same in-transition sentinel as Step 4.
                matchCorr = halfStep*(sameDirection&(prevTotalStep!=0));
                speedRatio = prevSpeed/(totalStep*invCurveScale*step+(1-active)*1e-30);
                clampedRatio = max(0, min(speedRatio, maxDerivVal));

                gonnaDo(phase) = (1-fdOf(phase))*totalStep;

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio))+matchCorr:min(1))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                forwardPos = select2(releasing,
                    // Attack branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    // Release branch
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        inverseDerivativeTopRelease(shape, clampedRatio)))+matchCorr:max(0):min(1-step);

                // The mirror branch: an attack whose inherited speed still
                // points UP is a turnaround. The matched anchor is the
                // negated EARLY-branch inverse — the segment of the attack
                // curve between phase 0 and the derivative peak, walked
                // toward the apex with speed decreasing to zero. (In the
                // inverse-derivative naming, which is inherited from the
                // release curve through the 1-x mirror, that early attack
                // branch is inverseDerivativeTOPAttack; BottomAttack is the
                // late branch that forward matching uses.) The anchor is
                // floored at the remaining budget so the landing can never
                // run past the arrival. On the flip sample itself
                // (direction just changed) the anchor is forced to the grid
                // depth -(arrival - att_samples)*attStep: the gate has just
                // certified that the matched depth reached it (to within
                // the one-sample crossing slack), and anchoring on the grid
                // is what makes the landing exact. Mid-turnaround
                // re-matches (retargets, clamp disturbances) use the
                // speed-matched depth, which by then can only be shallower
                // — landing early, never late.
                turning = attacking&(prevSpeed>0);
                mirrorMatched = max(
                    (0-inverseDerivativeTopAttack(shape,
                        min(max(0-speedRatio, 0), maxDerivVal)))+matchCorr,
                    0-corrBudget);
                mirrorPos = select2(sameDirection, 0-corrBudget, mirrorMatched);

                matchedPos = select2(turning, forwardPos, mirrorPos);

                // --- Step 6: advance phase ---
                //
                // anchor   = where we are on the curve NOW
                // newPhase = where we will be after this sample
                //
                // The 1-step clamp on the anchor makes the phase reach 1 (to
                // within rounding) on the final sample of a transition; the
                // landing test in Step 7 has half a step of slack to absorb
                // that rounding. Mirror anchors are negative and pass the
                // clamp untouched.
                anchor = select2(trusted, matchedPos, prevPhase:min(1-step));
                newPhaseRaw = (anchor+step)*active;

                // --- Step 7: curve increment + guarded landing ---
                //
                // The increment integrates the curve derivative over
                // [anchor, anchor+step] with 2-point Gauss-Legendre
                // quadrature on the cheap rational derivative (see v0.4
                // note B; the derivative of fdOf is the same expression for
                // both directions — the attack flip's two sign changes
                // cancel in the chain rule). derivM extends it to the
                // mirror: even position extension means an odd derivative,
                // signum(p)*fd'(|p|), which is LINEAR through the apex
                // (fd' ~ x near 0), so an interval straddling 0 integrates
                // cleanly — 2-point GL is exact through cubics.
                glA = 0.21132486540518713;
                glB = 0.78867513459481287;
                derivOf(p) = derivativeBaseRelease(shape,
                    select2(releasing, 1-p, p))*invCurveScale;
                derivM(p) = ma.signum(p)*derivOf(abs(p));
                delta = totalStep*step*0.5*(derivM(anchor+glA*step)+derivM(anchor+glB*step));

                // The value this loop wants to write. Stored as state so the
                // next sample can detect outside modification (see Step 4).
                expected = prev+delta;

                // Landing: when the phase completes, snap to the target —
                // but only across a residual the snap was designed for
                // (quadrature epsilon, float epsilon). guardSize is one
                // peak-speed sample. Under mid-transition retargets the
                // phase can complete with a real gap left (see v0.4 note
                // A); snapping across that is a teleport, so instead the
                // state resets WITHOUT snapping and the next sample starts
                // a fresh shaped transition over the leftover distance,
                // entered through velocity matching so it continues from
                // the current speed. The same refusal is what parks a HELD
                // release whose target dropped away: residual huge, no
                // snap, state resets, idle until the gate trips.
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
                // idle, and during a turnaround it equals lookaheadX, which
                // the rising mirror leg is above by construction. It cannot
                // conflict with the upper clamp because
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
//    Progress through the current transition. A release lives on [0, 1];
//    an attack on [-peakPhaseAttack, 1], where the negative leg is the
//    time-mirrored deceleration of a turnaround and 0 is the apex (zero
//    speed). Advances by `step` per sample; trusted between
//    discontinuities, re-derived from the measured speed (velocity
//    matching) across them.
//
//  turnaround:
//    The smooth reversal performed when an attack interrupts a release.
//    The attack curve is extended to negative phase by even reflection of
//    the fraction-done function: same positions, mirrored velocities. The
//    envelope rises decelerating along the mirror, stands still at the
//    apex, and falls through the regular attack — one C1 trajectory, no
//    corner. Sized by the `turnaround` slider (0 = v0.4 behavior, 1 = a
//    full extra attack of lookahead, enough for any matched depth).
//
//  arrival:
//    Countdown (samples) until the pending window minimum reaches
//    delayedX. Armed to lookaheadSamples when lookaheadX drops below the
//    envelope while nothing is pending; the first deadline binds, later
//    drops ride along as retargets.
//
//  held:
//    Wanting to attack while there is still time not to. The flip waits
//    until arrival <= att_samples + depth*(att_samples+1), so the landing
//    coincides with the arrival; meanwhile the in-flight release continues
//    on a frozen span, or the envelope idles.
//
//  totalStep:
//    The full amplitude span (in linear gain) of the current transition.
//    Maintained incrementally: held exactly while the target is static
//    (and during a HOLD), moved by exactly (lookaheadX - lookaheadX')
//    when it isn't, recomputed on entry as (lookaheadX - prev) — divided
//    by (1 - fd(grid depth)) for a turnaround entry, since the mirror leg
//    retraces ground — and never allowed to shrink past the span whose
//    curve can still carry the current speed (supportTotalStep, two-sided
//    on the attack).
//
//  step:
//    Phase increment per sample = 1 / (duration_in_samples).
//
//  fdOf(p):
//    Fraction of the transition completed at phase p, for both directions.
//    Position on the curve = transitionStart + totalStep*fdOf(p); for
//    mirror phases, fdOf(|p|). The per-sample delta integrates
//    signum(p)*fdOf'(|p|) over one step with 2-point Gauss-Legendre
//    quadrature; the ~1e-10 per-transition residual is absorbed by the
//    landing snap.
//
//  velocity matching:
//    Across a discontinuity (retarget, direction flip, outside
//    modification), find the phase on the (possibly new) curve whose
//    derivative equals the envelope's current speed, then continue from
//    there. A speed whose sign opposes the span maps to the mirror: minus
//    the early-branch inverse (TopAttack — the naming mirrors the release
//    curve), floored at the remaining budget. With
//    exact increments the measured speed is a mean over the last step,
//    hence the halfStep matchCorr.
//
//  gonnaMakeIt:
//    "If we re-anchor on the decelerating branch now, will the resulting
//    trajectory still reach the target?" Picks the accelerating vs
//    decelerating branch of the inverse derivative (forward matching
//    only; the mirror always decelerates into the apex).
//
//  landing:
//    When the phase completes, the transition state is fully reset, and
//    the output snaps to the target only if the residual is at most one
//    peak-speed sample (quadrature/float epsilon territory). Larger
//    residuals are NOT snapped (teleport); the leftover gets a fresh
//    shaped transition that velocity-matches into the current speed. The
//    phase-completion test carries halfStep of slack because the final
//    phase value is produced by accumulation and can round below 1. A
//    turnaround flipped on the grid lands on the exact sample the arrival
//    countdown reaches zero — the sample delayedX presents the transient.
//
//  lookaheadX:
//    The sliding minimum of the input over the next lookaheadSamples
//    (= att_samples + turnaround budget): the envelope's target.
//
//  AUC compensation:
//    (Omitted from this core file for clarity.) Adjusts the attack/release
//    duration so that changing the shape doesn't change the average level
//    of the envelope — sharper shapes get proportionally shorter durations.
