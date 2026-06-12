declare name "shapedSmoother_core";
declare version "0.14";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026 - 2026, Bart Brouns";

// ============================================================================
//  shapedSmoother — core algorithm, v0.14
//
//  CHANGES vs v0.4: THE TURNAROUND (v0.5), RESTRUCTURED FOR CPU (v0.6),
//  TURNAROUND FIXED AT 1 + SHARED INVERSION + SUPPORT RECIPROCAL (v0.7),
//  HEADING-ARMED COUNTDOWN — STAIRCASE-PROOF SCHEDULING (v0.8),
//  RELEASE BRANCH PICK — NO MORE PHANTOM-SPAN LUNGES (v0.9),
//  APEX CAPPED AT THE PRE-EPISODE MINIMUM (v0.10),
//  POSITION-BUDGET DIVISION OFF THE RECURRENCE PATH (v0.11),
//  TURNAROUND RETIRED — THE BRAKE (v0.12),
//  BRAKE RE-CLOCKED OFF THE FUTURE TARGET (v0.13),
//  THE BRAKE FINISHED — SETTLE, COUNTDOWN, REACH RULE, AND A DIAL (v0.14).
//
//  v0.12 retires the turnaround. Its guarantee had narrowed with every
//  correctness fix (over-capacity speeds and ceiling-bound rises both
//  velocity-dump), its corrected form had already given back most of
//  the level retention that motivated it, and its scheduling machinery
//  — countdown, heading-armed gate, mirror branch, apex ceiling,
//  position budget — was the parent of every late-corner bug class
//  (v0.8–v0.10). The brake (Step 1, after v0.3.6) replaces all of it:
//  a second sliding minimum sees one brake window further than the
//  attack window, and while a deeper minimum is visible there the
//  release increment ramps linearly to zero, so the attack always
//  enters from a standstill through ordinary velocity matching — C1
//  at any inherited speed, punctual by the plain v0.4 contract
//  (attacks start the moment the attack window binds and take
//  att_samples+1), and structurally incapable of sitting above a
//  level still in flight. Latency is unchanged: att + (att+1), now
//  spent on warning instead of reversal. Measured on the 60 s noisy
//  staircase panels: ZERO deep corners at both; the average envelope
//  sits 0.026–0.039 HIGHER than v0.11, because the release now
//  targets the attack-window minimum instead of the doubled-window
//  minimum, so near-term climbs survive; faust compile 3 s (from
//  25 s); CPU 0.81x of v0.11 under -double -mcd 0 -Ofast. Paragraphs
//  below this one describe the retired versions and are kept as
//  history.
//
//  v0.13 re-clocks the brake. v0.12 started the ramp when the
//  envelope lifted above the brake horizon, so a release born shortly
//  before a bind could reach the flip with the ramp only partially
//  built and carry up to a third of peak speed into it — a visible
//  slope corner. The ramp is now clocked by the gap between the
//  FUTURE attack target (the same sliding minimum read one brake
//  window earlier in its pipeline, costing nothing) and the present
//  one. That gap opens exactly brake_samples before every bind and
//  closes at it, so the ramp saturates exactly when the attack window
//  binds, for every release however young. Consequence: a release
//  born with the gap already open starts already-braked, so
//  recoveries shorter than the brake window are suppressed instead of
//  blipped — the policy a release hold applies. The second
//  sliding-minimum tree is gone (the gap needs no separate horizon
//  minimum). Measured: the scope's release-into-bind corner (36% of
//  peak speed at the flip) is eliminated for standstill births and
//  for fresh deeper minima, whose ramp timing is exact; equal-depth
//  successions caught mid-flight stay partially braked — worst flip
//  18% of peak on the staircase loop, half of v0.12 — since exact
//  timing there needs per-minimum age, the countdown machinery this
//  design retired. Both noise panels: zero deep corners. Average
//  envelope on dense noisy material pays for the always-armed brake:
//  -0.342/-0.287 vs v0.12's -0.230/-0.224 — sub-brake recoveries are
//  suppressed by design (a release hold does the same). CPU 83.4 ns
//  per sample (1.08x of v0.12, generic -double -mcd 0 -Ofast).
//
//  v0.14 finishes the brake, prices it, and hands over the dial.
//  The net of one long bench session, in four parts:
//
//  SETTLE ONTO LEVELS, RIDE SLOPES. v0.13 froze every braked
//  release, including the harmless ones below the coming floor.
//  Those now retarget to the floor itself and run an ordinary full
//  shaped release INTO the level the envelope will be asked to hold
//  — brisk through the middle, decelerating only in the shape's own
//  tail, landing and parking. (A velocity-cap variant settling
//  exponentially over one brake window was tried first and replaced:
//  the asymptote never lands, and dense material crawled at a few
//  percent of curve speed.) The retarget fires only while the floor
//  is STEADY — a reigning deepest sample holds attackMinAhead
//  bit-constant, a slope's minimum drifts every sample — so smooth
//  descents are ridden at the attack window, not slope*(att+brake)
//  under the signal (0.4 of level at the worst smooth-descent
//  motif); approaching from below can never end above the floor, so
//  every brake guarantee survives.
//
//  RISE WHILE YOU CAN — THE RAMP IS A COUNTDOWN. The brake had
//  grown into a blanket ban on any rise it would have to give back
//  within the horizon. But rise-and-duck is C1-safe whenever the
//  rise reaches ZERO SPEED by the bind, and the bind is knowable:
//  when a gap opens by attackMinAhead DROPPING, the bind is exactly
//  brake_samples later, so the ramp doubles as a countdown — time to
//  bind = (1-brakeRamp)*brake_samples — and any later, deeper drop
//  only pushes the true bind further out (gapFresh). An above-floor
//  release whose remaining flight fits a trusted countdown runs
//  FREE: it lands at the present target and the attack departs from
//  a standstill AT the target — the full bump. One that cannot land
//  in time fades linearly to zero exactly at the bind — a partial
//  bump, still kink-free. Fading BELOW-floor releases too was tried
//  and reverted: perpetual gaps on dense material keep the ramp
//  saturated and freeze releases that threaten nothing (the v0.13
//  level loss, reimported). Only countdown-invalid standstill
//  births remain held.
//
//  THE REACH RULE RESTORED. v0.9's release branch pick
//  (strand-versus-overshoot) existed to duck the upward-retarget
//  lunge of the turnaround era; its overshoot estimate overstates
//  late in flight, so dense upward retargets kept re-anchoring
//  rising releases onto the decelerating tail — the envelope
//  rounded off and parked far below its target, 0.03-0.04 of mean
//  envelope in BOTH dial positions. v0.4's reach rule is back, and
//  the lunge is fixed at its root: while braking, an upward release
//  retarget re-derives the span entry-style (Step 3), restoring
//  span/position consistency by construction; outside braking the
//  surviving lunges measure at v0.4 parity (max ~0.004, zero deep
//  corners) — the cost the reference behavior always paid, caught
//  by the ordinary attack-back.
//
//  THE DIAL. Smoothness and level retention turned out to be ONE
//  trade with two ends — the brake's corner-free binds against
//  v0.4's full-speed bumps — so brakeEnable (below) selects the end
//  instead of the algorithm choosing for you. Measured on the full
//  battery (two 60 s noisy staircase panels, the blocky and noisy
//  user scenes, the square loop), zero deep corners everywhere in
//  both modes. ON: panel mean envelope -0.2075/-0.2097 (v0.13:
//  -0.342/-0.287), worst flip 0.0004-0.0014 at 4-26 per minute, the
//  square loop at zero flips; the residual flips are the
//  masked-succession class, bounded by the per-minimum-age
//  limitation. OFF: bit-equal to the reach rule on the v0.4
//  baseline — mean envelope within 0.002 of v0.4 on every scene,
//  with v0.4's flip statistics (600-1300 per minute at 40-77% of
//  peak: the corners the brake exists to remove) — and the brake
//  window of latency collapses with it. The open reconciliation —
//  releases running free until the latest moment a curve-tail
//  deceleration can still reach zero by the bind — needs
//  capped-target machinery and is future work.
//
//  An attack arriving while a release was still in flight hit the
//  clampedRatio floor: the opposite-sign speed matched to phase 0, the
//  envelope stopped dead in one sample and re-accelerated from rest — an
//  instant change of direction. v0.5 bought extra lookahead and spent it
//  on turning around smoothly; v0.6 is the same feature with the two
//  audio-rate transcendentals it introduced designed back out. Per sample
//  the loop now costs one rational evaluation, one sqrt, a counter and
//  abs/signum over v0.4 (measured numbers at the bottom of this header).
//
//  v0.7 fixes the turnaround budget at 1 (the slider is gone — see the
//  note above extraSamples; bit-exact, since the slider's default was
//  already 1) and trims the audio loop twice. First, the forward and
//  mirror inverse-derivative inversions are consumed mutually
//  exclusively (matchedPos selects on `turning`), so the ratio is
//  direction-folded BEFORE the inversion — a sign flip plus clamp,
//  exact in float; see Step 5 — and one sqrt and one rational
//  denominator serve both paths at bit-identical output. Second, the
//  support floor multiplies by slider-rate reciprocals instead of
//  dividing (ulp-level, Step 3) — that division sat on the
//  recurrence's latency chain and its removal is most of the measured
//  speedup.
//
//  v0.8 fixes the turnaround sometimes arriving late. The countdown
//  used to arm only on drops below the CURRENT envelope; a release
//  interrupted by a descending staircase would climb — held span,
//  rise leg — above intermediate input levels that had entered the
//  window above the envelope, and those levels reached the output tap
//  with only the hard clamp to duck them: a corner, scheduled for
//  nobody. The countdown now arms against the envelope's HEADING (the
//  live release's target; Step 0), so every level the envelope could
//  climb above carries its own deadline, and the hold then freezes
//  the span below everything still in flight. Exact per-horizon
//  scheduling — the ramp-plus-sliding-min of linear-attack limiters —
//  is not available here: a fixed-duration shaped attack decides its
//  span at the flip, so a sample's pass-under time cannot be priced
//  at window entry. The heading rule closes the reachable class
//  instead; what remains is the one-sample crossing slack and — the
//  dominant residue — the rise leg's apex, which is placed by the
//  INHERITED SPEED alone and can still overshoot levels in flight
//  (capping it at the frozen heading would close that class too, at
//  the price of a bounded velocity mismatch and one more state word;
//  landed in v0.10). Measured on 60 s of noisy staircases
//  (noise 0.152, blockscale 8.21, att 4.5 ms, rel 14.4 ms,
//  shapes 0.5/0): clamp corners at transient arrivals 106 -> 28,
//  corners deeper than 0.01 of full scale 81 -> 15, worst corner
//  duration 8 -> 2 samples; cost 0.011 of extra average attenuation
//  and 1.9% more fully-ducked time.
//
//  v0.9 removes the worst of what remained. Velocity matching's
//  branch pick preferred the anchor that REACHES the target in both
//  directions; on the release side, after an upward retarget near
//  completion, that pick re-anchors a grown span at its start and
//  lunges past the target by the already-traversed distance, riding
//  the input clamp and parking far above every level still in flight
//  — the setup for the deepest late-turnaround corners. The release
//  now reaches only when the decelerating anchor's strand exceeds the
//  accelerating anchor's overshoot (see Step 5); the rule is
//  unchanged at fresh entries. Measured on the same 60 s of noisy
//  staircases, att 7.41 ms / shapes 0.341 / 0.5: worst corner depth
//  0.466 -> 0.167, corners 38 -> 24, all residuals one sample; the
//  average envelope sits ~0.02 lower because the lunge had been an
//  illegitimately fast release. What remains is the speed-placed
//  apex class described above.
//
//  v0.10 lands the cap. One new recursion word latches the window
//  minimum from the last sample before a hold or attack episode
//  begins; every input sample arriving before the episode's landing
//  was inside the window at that freeze, so it sits at or above the
//  latch, and a rise apex capped there can no longer overshoot
//  anything still in flight. The cap bounds the anchor depth using
//  fd(q) <= q*maxDerivativeBase*invCurveScale — conservative, so the
//  trade is always a small velocity mismatch, never a position
//  corner — and a capped rise only lands earlier, leaving the
//  schedule's punctuality untouched. The hold itself also ends when
//  the held release — possibly riding a floor-inflated span past its
//  pre-drop target — reaches the same latch, so the climb can never
//  pass a level still in flight either. Measured on the 60 s noisy
//  staircases at both panels (att 4.51 ms, shapes 0.5/0; and att
//  7.41 ms, shapes 0.341/0.5): corners deeper than 0.01 of full
//  scale: ZERO, at both — every residual is the one-sample crossing
//  slack — with rise counts and average attenuation within 1% of
//  v0.9.
//
//  1. MIRRORED ATTACK CURVE / NEGATIVE PHASE. The attack's fraction-done
//     function is extended to negative phase by even reflection,
//     fd(p) := fd(|p|), so position(p) = apex + totalStep*fd(|p|), and the
//     chain rule flips the velocity's sign across zero: the derivative of
//     the extension is signum(p)*fd'(|p|). Advancing the phase from -q
//     through 0 traces a rise that decelerates to a standstill at the
//     apex (phase 0). The velocity is continuous everywhere, including
//     the apex (zero from both sides; fd is locally parabolic there, so
//     the 2-point quadrature crosses it without noticing). Cost: abs() on
//     the curve evaluations, signum() on the derivative.
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
//     Velocity matching would place the current speed at mirror depth d
//     with derivative(d) = speed; the rise from -d plus the attack proper
//     costs att_samples + d*(att_samples+1) samples, so the flip waits
//     until the countdown has shrunk to exactly that. Rather than
//     inverting the derivative (a square root and several divisions per
//     sample in the first draft), the gate compares ON THE DERIVATIVE
//     AXIS, which is monotone over the early branch:
//
//         hold while  budget > peak,
//             or while  derivative(budget)*|span| > max(speed, 0)
//
//     where budget = (arrival-att_samples)*attStep. Depths beyond the
//     derivative peak cannot be velocity matched, hence the first term; a
//     speed above the curve's capacity flips once the budget shrinks to
//     the peak and carries the bounded mismatch. The comparison flips on
//     the same sample as the inverted form — exactly — and costs one
//     rational evaluation.
//     Downward or zero speed degenerates to d = 0: the flip waits for
//     arrival == att_samples, v0.4's start time to the sample. The
//     gate only ever delays ENTRY into the attack (prevTotalStep >= 0):
//     once flipped, the schedule is committed and an attack in flight is
//     never re-gated — the rise sits exactly on the gate's own equality,
//     so re-testing it would reduce to rounding noise. (Up to v0.6 a
//     `turnaround` slider scaled the budget, and 0 made the feature
//     inert — with exactBranchPick = 1 that was a bit-for-bit A/B
//     against v0.4. v0.7 fixes the budget at 1, so that escape hatch is
//     gone; exactBranchPick remains an independent exact-vs-fitted
//     switch for the landing projection, see Step 5. The default
//     approximated branch pick can resolve razor-edge matching
//     recoveries differently.)
//
//     During the hold the in-flight release simply continues — its span
//     frozen against the (already dropped) target, its phase trusted —
//     and if it completes it parks: the landing guard refuses the snap
//     and the state resets to idle until the gate trips. A release too
//     fast for the attack to match therefore decelerates to a standstill
//     on its OWN curve, idles, and the attack launches on schedule:
//     corner-free by a different route.
//
//  4. TWO LEGS, APEX HANDOFF. The turnaround is two transitions, not one.
//     The RISE leg anchors at the grid depth -(arrival-att_samples)
//     *attStep with the plain entry span (lookaheadX - prev) — the same
//     span the gate matched against, so the velocity at the flip is exact
//     to within the one-sample crossing slack. The rise needs no position
//     exactness at all: reaching phase 0 is a leg COMPLETION — the state
//     resets like any completed transition, except it can never snap (the
//     target is a full attack away). The next sample re-enters the attack
//     through ordinary velocity matching, with the apex's near-zero speed
//     anchoring it at phase ~0 and a fresh span equal to the TRUE
//     remaining distance. Landing position and landing time are exact by
//     construction — rise (q/attStep samples) plus descent
//     (att_samples+1) meets the countdown — with no curve evaluation.
//     (The first draft instead stretched a single span by 1/(1-fd(q)),
//     which cost a log and an atan every sample and left a bounded slope
//     error at the flip because the gate had matched the unstretched
//     span. Both are gone.) The only residue of the handoff is an
//     acceleration step at the apex proportional to the difference of
//     the two spans — at most about half the trajectory's own peak
//     acceleration in the worst case, usually orders below it, and a
//     class gentler than the velocity step it replaces.
//
//  5. SENTINEL MIGRATION. Phase 0 is now a legitimate value DURING a
//     turnaround (the rise walks up to it), so "no transition in
//     progress" moves from prevPhase == 0 to prevTotalStep == 0 — a value
//     the completion reset already establishes and which no live
//     transition holds. The trusted gate and the matchCorr gate use the
//     new sentinel.
//
//  6. TWO-SIDED SUPPORT FLOOR (attack). During the rise the speed has the
//     opposite sign to the span, so the attack-side floor uses
//     -abs(prevSpeed): a span born in a turnaround can always carry the
//     speed it inherits. The release side keeps the signed prevSpeed —
//     entry from an attack must stay inert there.
//
//  COSTS. Latency rises from att_samples to
//  att_samples + (att_samples+1): the turnaround budget is one full
//  extra attack, i.e. doubled lookahead. Per sample, measured against
//  v0.4 in the generated C++ (48 kHz, -O3 -ffast-math, block 256,
//  random steps): v0.4 1.00x, first turnaround draft 2.4x, v0.6 1.41x
//  — or 1.86x with exactBranchPick = 1. v0.6 measured against THIS
//  version, same flags, isolated smoother (input-driven process, no
//  test-signal cost), Xeon 2.8 GHz: v0.6 = 1.05–1.10x of v0.7
//  across repeated runs, at the defaults and at att = 2 ms /
//  shapes = 0.9 alike. v0.8's arming change costs one max and one
//  select2 per sample and no state (corner statistics in the v0.8
//  paragraph above); v0.10 adds one recursion word and one division
//  on the mirror path (the position budget, Step 5), and v0.11 feeds
//  that division from state and feed-forward signals only, so it
//  issues at sample start and overlaps the recurrence instead of
//  lengthening it — on narrow cores the difference is real. Nearly all of the v0.7
//  speedup is the support-floor reciprocal (Step 3); the shared
//  inversion is
//  time-neutral on a wide core (the discarded twin ran in parallel on
//  the divider port) but halves the loop's sqrt count at bit-identical
//  output. Audio-rate transcendentals per sample: v0.4 = 1 log,
//  1 atan, 2 sqrt; this version's default = 0 log, 0 atan, 1 sqrt —
//  the forward and mirror inversions share one evaluation (Step 5),
//  and the branch pick reads slider-rate-fitted cubics instead of
//  evaluating the curve. (exactBranchPick = 1 brings the log+atan pair
//  back: it is v0.4's own gonnaMakeIt projection, kept verbatim; its
//  per-shape constants are hoisted, bit-identically, through the
//  direction select.) What the turnaround itself adds is rational and
//  cheap: the gate's cross-multiplied comparison, the arrival counter,
//  signum/abs in the quadrature, and one more sliding-min level from
//  the doubled window. The lookahead clamp keeps the generated ring
//  memory at v0.4's footprint.
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
relMs = SmootherGroup(hslider("[3]rel[scale:log]", 0.05*1000, 1, 5000, 0.1));
releaseShapeSlider = SmootherGroup(hslider("[4]release shape", 0, 0, 1, 0.001));

// Derived parameters
maxHold = 0.05;
maxSR = 48000;
maxHoldSamples = maxHold*maxSR;
maxLookaheadSamples = 2*maxHoldSamples+1;

// Compile-time switch. 1 makes the landing projection inside velocity
// matching evaluate the exact curve — v0.4's own projection, verbatim
// (the bit-for-bit v0.4 A/B it used to enable required turnaround = 0,
// which v0.7 removed) — at the price of an audio-rate
// log+atan pair, the single most expensive thing in the loop (see
// COSTS). The default evaluates per-direction cubics fitted to the
// curve's deceleration tail at slider rate instead; the pick they
// produce can differ from the exact one only when the projected landing
// sits within the fit error of the target (2e-4 of the span for the
// attack, 1.3e-2 for the release), where either pick is absorbed by the
// landing machinery. Faust folds the constant, so the unused path costs
// nothing.
exactBranchPick = 0;

att = attMs/1000;
rel = relMs/1000;
att_samples = att*ma.SR:max(1);
rel_samples = rel*ma.SR:max(1);

// Brake budget (v0.12). The brake horizon extends the lookahead by
// brake_samples beyond the attack window, purely so an upcoming attack
// is visible early enough for the release to glide to a stop before
// the attack window binds. One full extra attack keeps the total
// latency identical to the turnaround versions (v0.5–v0.11); sizing it
// larger (e.g. min(maxBrakeSamples, max(att, rel) in samples)) brakes
// more gently at more latency — the deceleration magnitude is
// speed/brake. lookaheadSamples is the single source of truth for the
// delay-line length and the latency.
// THE DIAL (v0.14). brakeEnable = 1 is the brake: every release
// arrives at every bind at ~zero speed (worst measured flip
// 0.0004-0.0014 across the test battery, square loop zero), at a
// measured cost of 0.017-0.026 of average envelope and partial
// recovery bumps. brakeEnable = 0 is the v0.4 trajectory family:
// full-speed rises, full bumps, mean envelope within 0.002 of v0.4
// on every scene — and a velocity corner at every release-to-
// attack flip (600-1300 per minute on dense material, up to 77% of
// peak speed). Punctuality (zero deep corners) holds in BOTH
// modes; the choice is purely corner-versus-level. With
// brakeEnable = 0 every brake mechanism strips to nothing and the
// extra brake window of latency collapses with it:
// lookaheadSamples falls back to att_samples and attackMinAhead
// coincides with lookaheadX, so braking is constant-false and the
// compiler removes the rest.
brakeEnable = 1;
brake_samples = (int(att_samples)+1)*brakeEnable;
brakeStep = 1/max(brake_samples, 1);
lookaheadSamples = min(int(att_samples)+brake_samples, int(maxLookaheadSamples));

process = testSignal<:(_@lookaheadSamples, shapedSmoother(_));

shapedSmoother(x) = lookaheadX, attackMinAhead, delayedX:env~(_, _, _, _):(_, _, _, !)
    with {
        // Future attack window: the att_samples+1 sliding minimum, not
        // yet delayed — it equals lookaheadX brake_samples in the
        // future, for free. A sample entering it hands off to
        // lookaheadX exactly brake_samples later (brake_samples =
        // att_samples+1 = the window length, so the two windows tile
        // the lookahead without a gap; a differently sized brake would
        // need the brake-only band's own minimum here instead).
        attackMinAhead = x:ba.slidingMin(int(att_samples)+1, 1+maxLookaheadSamples);
        lookaheadX = attackMinAhead@brake_samples;
        // THE BRAKE CLOCK (v0.13). braking is true exactly while the
        // future target is below the present one — from the moment a
        // deeper minimum becomes visible until the moment the attack
        // window binds on it, a span of exactly brake_samples — so the
        // ramp saturates exactly at every bind and any release in
        // flight is at a standstill when its attack begins, however
        // recently the release was born. A release born with the gap
        // already open starts already-braked: recoveries shorter than
        // the brake window are suppressed instead of blipped, the same
        // policy a release hold applies. (v0.12 clocked the ramp from
        // the envelope's own liftoff, which could leave it only
        // partially built at the bind — a release could carry a third
        // of peak speed into the flip.)
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

        releasePeakPhase = peakPhaseRelease(releaseShape);

        // Tail fit (for the default approximated branch pick, Step 5):
        // per direction, two cubics interpolating the exact fraction-done
        // curve over the deceleration tail [peakPhase, 1] — the only
        // region the landing projection ever evaluates — split at
        // u = 0.15 in tail coordinates. Fitted HERE, at slider rate,
        // where the curve's log+atan cost nothing; the loop then sees
        // only a Horner evaluation. Max fit error across the full shape
        // range: 2e-4 for the attack tail, 1.3e-2 for the release tail
        // (whose tail spans nearly the whole curve at sharp shapes).
        tailSplit = 0.15;
        // monomial coefficients of the cubic through f at t = 0, 1/3, 2/3, 1
        tc0(f0, f1, f2, f3) = f0;
        tc1(f0, f1, f2, f3) = 9*f1-5.5*f0-4.5*f2+f3;
        tc2(f0, f1, f2, f3) = 9*f0-22.5*f1+18*f2-4.5*f3;
        tc3(f0, f1, f2, f3) = 13.5*f1-4.5*f0-13.5*f2+4.5*f3;
        attackTailFd(u) = cheapCurveAttack(attackShape,
            attackPeakPhase+u*(1-attackPeakPhase));
        releaseTailFd(u) = cheapCurveRelease(releaseShape,
            releasePeakPhase+u*(1-releasePeakPhase));
        attLoF(k) = attackTailFd(k*tailSplit/3);
        attHiF(k) = attackTailFd(tailSplit+k*(1-tailSplit)/3);
        relLoF(k) = releaseTailFd(k*tailSplit/3);
        relHiF(k) = releaseTailFd(tailSplit+k*(1-tailSplit)/3);

        // Fencepost: a transition occupies the CLOSED interval of N+1 output
        // samples [entry .. entry+N], so the phase advances in N+1 steps of
        // 1/(N+1). With 1/N the attack landed one sample early: lookaheadX
        // drops at m (the new minimum enters the window immediately) while
        // delayedX delivers the transient at m+480, and N increments
        // starting ON the entry sample complete at m+479. With 1/(N+1) the
        // attack's landing snap fires on exactly the sample where delayedX
        // presents the transient, and a release elapses exactly
        // rel_samples from onset to target. (The same arithmetic is what
        // makes the rise + descent of a turnaround meet the arrival
        // countdown of Step 0 on the sample.)
        attStep = 1/((att*ma.SR):max(1)+1);
        relStep = 1/((rel*ma.SR):max(1)+1);

        // Reciprocal of the per-direction support-floor denominator
        // (maxDerivativeBase*invCurveScale*step), hoisted to slider rate
        // so Step 3 multiplies instead of dividing (v0.7): that division
        // sat on the recurrence's latency chain, and it is most of
        // v0.7's measured speedup. Within ~1.5 ulp of the divided form
        // per sample; like any non-bit-exact change inside the loop it
        // can resolve a razor-edge branch pick differently and shift one
        // leg's trajectory transiently before the landing machinery
        // re-converges it (measured on 30 s of random steps: <= ~1e-4 of
        // full scale at the defaults, <= ~3e-3 at att = 2 ms /
        // shapes = 0.9) — the same class of difference the default
        // branch pick already carries vs exactBranchPick = 1.
        attackSupportRecip = 1/(attackMaxDerivBase*attackInvCurveScale*attStep);
        releaseSupportRecip = 1/(releaseMaxDerivBase*releaseInvCurveScale*relStep);

        // =======================================================================
        //  env: the recursive core
        //
        //  State carried between samples:
        //    prev          — current envelope value
        //    prevPhase     — where we are on the curve. [0, 1] for a release
        //                    or a plain attack; [-attackPeakPhase, 0] for the
        //                    rise leg of a turnaround, whose negative phase
        //                    addresses the mirrored attack curve.
        //    prevTotalStep — total span of the current transition. Also the
        //                    idle sentinel: == 0 exactly iff no transition is
        //                    in progress (the completion reset establishes it,
        //                    and phase 0 can no longer serve — the rise leg
        //                    ends there).
        //    prevExpected  — the value THIS loop computed last sample.
        //                    prev != prevExpected means something outside the
        //                    loop (the output clamps, or an external
        //                    algorithm) modified the envelope, so the stored
        //                    phase can no longer be trusted and we must
        //                    velocity-match.
        //    prevArrival   — samples until the pending window minimum reaches
        //                    delayedX; <= 0 when nothing is pending.
        // =======================================================================
        env(prev, prevPhase, prevTotalStep, prevExpected, lookaheadX, attackMinAhead, delayedX) = result, newPhase, totalStepOut, expected
            with {
                prevSpeed = prev-prev';


                // --- Step 1: direction, gated by the turnaround schedule ---
                //
                // Wanting to attack no longer means attacking NOW. Velocity
                // matching would place the current speed at the mirror
                // depth d whose derivative equals it; the rise from -d plus
                // the attack proper costs att_samples + d*(att_samples+1),
                // so the flip waits until the countdown equals that. The
                // budget the schedule can still afford is
                //     mirrorBudget = (arrival - att_samples)*attStep,
                // and instead of inverting the derivative to find d, the
                // gate compares on the derivative axis, which is monotone
                // over the early branch (budget <= matched depth iff
                // derivative(budget) <= speed in span units). The
                // comparison flips on the same sample as the inverted form,
                // and it is cross-multiplied against the curve denominator
                // (always positive: c*(x-1)^2+(1-c) >= 1-c > 0), so the
                // loop carries a handful of multiplies and no division for
                // it. Depths beyond the
                // derivative peak cannot be velocity matched, so the hold
                // also persists while the budget still exceeds the peak: a
                // speed above the curve's capacity flips only once the
                // budget has shrunk to the peak, carrying the bounded
                // mismatch (in practice such releases usually park first —
                // see below).
                // Downward or zero speed degenerates to d = 0: the flip
                // waits for arrival == att_samples, v0.4's start time to
                // the sample.
                //
                // The gate is an ENTRY gate: it may only delay the flip
                // while no attack span is in flight (prevTotalStep >= 0 —
                // idle, or a release to hold). Once the flip has happened
                // the schedule is committed, and the gate is disarmed: the
                // rise traverses the very equality the gate tests — budget
                // and matched depth shrink in lockstep by construction —
                // so re-evaluating it there is a coin toss on rounding
                // noise, and one wrong sample wipes the transition state
                // and freezes the envelope until the budget drains (found
                // the hard way, in float32, at sharp shapes). The same
                // term protects a mid-attack retarget: a deeper minimum
                // arming a fresh countdown must not hold the attack it
                // lands in.
                //
                // During the hold the in-flight release (if any) continues:
                // span frozen (Step 3), phase trusted (Step 4). If it
                // completes first, the landing guard refuses to snap across
                // the dropped target and the state parks at idle until the
                // gate trips.
                attacking = (prev>lookaheadX)&(att>0);

                // The clock: a deeper-or-equal minimum is in the
                // brake-only region exactly when the future target sits
                // below the present one, and it binds exactly
                // brake_samples after appearing there.
                braking = attackMinAhead<lookaheadX;

                // COUNTDOWN VALIDITY (v0.14). The ramp doubles as a
                // bind countdown — time to bind = (1-brakeRamp)*
                // brake_samples — but only when the current gap's
                // deepest minimum announced itself by DROPPING into
                // the future window (then its bind is exactly
                // brake_samples after the drop, and any later, deeper
                // drop only pushes the true bind further out, so the
                // countdown stays safe-or-conservative). A gap opened
                // by the TARGET rising has a stale clock: the pending
                // minimum was already mid-pipeline, its bind is
                // sooner than the ramp claims, and nothing scheduled
                // on it can be trusted.
                gapFresh = (max(attackMinAhead<attackMinAhead')*braking)~_;

                // RELEASE RETARGET (v0.14). While braking, a release
                // BELOW the future floor no longer aims at lookaheadX
                // with a velocity cap — it aims at the floor itself.
                // The envelope then runs an ordinary full shaped
                // release INTO the level it will be asked to hold:
                // brisk through the middle, decelerating only in the
                // release shape's own tail, landing (snap, park)
                // instead of settling asymptotically. A velocity cap tried first was
                // an exponential law with the brake window as time
                // constant — C1, but decelerating everywhere and never
                // arriving; on dense material it crawled at a few
                // percent of curve speed with most of the headroom
                // still unused. The retarget keeps every guarantee (an
                // approach from below can never end up above the
                // floor, whatever the durations) and reuses the whole
                // existing transition machinery: target changes are
                // ordinary retargets, absorbed by velocity matching.
                // Above the floor nothing changes: the target stays
                // lookaheadX and the time fade glides the release to a
                // standstill. relTarget equals lookaheadX in every
                // attacking state (attacking means prev > lookaheadX
                // >= attackMinAhead), so attack spans are untouched.
                // ...but only onto LEVELS, not slopes (v0.14). A
                // genuine floor is one deepest sample reigning over
                // the future window, so attackMinAhead is bit-constant
                // while it reigns; on a gently descending input the
                // window minimum is the trailing edge and drifts every
                // sample. Without this gate, braking is true on ANY
                // descent (a later window's min is always lower on a
                // downslope) and the retarget pins the release to the
                // input's value one full horizon ahead — the doubled-
                // window tracking gap of v0.5–v0.11 returning through
                // the brake, slope*(att+brake) of level lost on every
                // smooth descent. Steady floors get the shaped settle
                // and its zero-kink guarantee; drifting slopes are
                // ridden at the attack window like v0.12, where the
                // envelope meets the moving target in its own shape
                // tail and tracks it through standstill micro-landings.
                // The one-sample detection latency of a fresh floor is
                // absorbed by velocity matching at the retarget.
                floorSteady = (attackMinAhead==attackMinAhead');
                relTarget = select2(braking&(prev<attackMinAhead)&floorSteady,
                    lookaheadX, attackMinAhead);

                // BIRTH HOLD (v0.13). The gap can also open by the
                // TARGET rising (a recovery window shorter than the
                // brake: the coming minimum equals the floor the
                // envelope already sits on, so attackMinAhead never
                // dropped). The bind is then less than a brake away and
                // no time ramp can saturate for it — so a release that
                // would be BORN into that situation simply isn't: a
                // true standstill (prevTotalStep == 0, exact — the
                // completion reset writes literal zero) stays parked
                // while a duck to at-or-below the current level is
                // pending. Velocity 0 -> 0: no kink is possible.
                // Recoveries shorter than the brake window are
                // suppressed, the policy a release hold applies.
                // Anything already in flight is untouched — the gradual
                // ramp below covers it, and a multiplier step at a
                // retarget is absorbed by velocity matching, which
                // re-anchors at the measured speed.
                freshHold = braking&(attackMinAhead<=prev)&(prevTotalStep==0)
                    &(1-gapFresh);
                releasing = (prev<relTarget)&(rel>0)&(1-freshHold);

                // THE BRAKE's time fade (v0.12, re-clocked v0.13).
                // The ramp integrates braking; it must integrate HERE,
                // where braking arrives as a plain input — written at
                // shapedSmoother level the recursion's feedback binds
                // to the wrong free input, x leaks in through
                // braking's definition, and the multiplier silently
                // leaves [0, 1]. The gap opens exactly brake_samples
                // before a fresh deeper minimum binds, so the fade
                // reaches zero exactly at the bind: an in-flight
                // release that cannot land in time (fading, Step 7)
                // glides to a standstill and its attack departs from
                // ~zero speed. Releases below the floor are handled by
                // the retarget above, not by the fade.
                //
                // Residual (bounded, documented): one scalar ramp and
                // one future minimum cannot time STAGGERED threats
                // individually. A deeper minimum arriving mid-settle
                // flips the regime from shaped release to saturated
                // fade — a deceleration step bounded by the current
                // curve speed — and an equal-depth succession caught
                // mid-flight keeps a partial fade. Exact per-threat
                // timing needs per-minimum age: the countdown
                // machinery this design retired.
                brakeRamp = ((_+brakeStep):min(1)*braking)~_;
                active = attacking|releasing;

                // --- Step 2: Select shape parameters for current direction ---
                shape = select2(releasing, attackShape, releaseShape);
                invCurveScale = attackInvCurveScale+releasing*(releaseInvCurveScale-attackInvCurveScale);
                zeroVal = select2(releasing, attackZero, releaseZero);
                maxDerivVal = select2(releasing, attackMaxDerivBase, releaseMaxDerivBase);
                step = select2(releasing, attStep, relStep);
                durSamples = select2(releasing,
                    (att*ma.SR:max(1))+1, (rel*ma.SR:max(1))+1);
                halfStep = 0.5*step;

                // --- Curve helpers ---
                // cbAt(p): normalized curve value (attack flip applied inside)
                // fdOf(p): fraction of the transition completed at phase p,
                //          0 at the start, 1 at the target, for BOTH directions.
                //          Position on curve = transitionStart + totalStep*fdOf(p).
                //          For negative (rise leg) phases the position uses
                //          fdOf(|p|); only the derivative needs the sign
                //          (see derivM in Step 7).
                // cheapCurveBase's per-shape constants are hoisted through
                // the direction select so the loop carries no sqrt or
                // reciprocal for them; the arithmetic order is preserved,
                // so the result is bit-identical to the unhoisted form.
                // (Alive only when exactBranchPick = 1 — otherwise this
                // whole chain is pruned.)
                kVal = select2(releasing,
                    sqrt((1/attackShape)-1), sqrt((1/releaseShape)-1));
                twoC = select2(releasing, 2*attackShape, 2*releaseShape);
                oneMinusC = select2(releasing, 1-attackShape, 1-releaseShape);
                cbAt(p) = ((log(shape*xx*xx+oneMinusC)
                    +2*kVal*atan(xx/kVal)-2*xx)/twoC-zeroVal)*invCurveScale
                with { xx = select2(releasing, 1-p, p); };
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
                // distance. The rise leg of a turnaround enters with this
                // plain span too: it needs no position exactness (the
                // descent re-spans at the apex, Step 7), and the gate
                // matched the speed against exactly this span, so the flip
                // is velocity-exact to within the crossing slack.
                //
                // supportTotalStep is the smallest span whose curve can
                // carry the current speed. The attack side uses
                // -abs(prevSpeed): during the rise the speed has the
                // opposite sign to the span, and a turnaround must be born
                // with a span that can carry the speed it inherits. The
                // release side keeps the signed prevSpeed so that entry
                // from an attack stays inert, as in v0.4. On the floor the
                // envelope decelerates along the curve at the curve's
                // steepest rate; the deliberate overshoot this allows is
                // caught by the output clamps in Step 7.
                incrementalTotalStep = prevTotalStep+(relTarget-relTarget');
                entryTotalStep = relTarget-prev;
                supportTotalStep = select2(releasing, 0-abs(prevSpeed), prevSpeed)
                    *select2(releasing, attackSupportRecip, releaseSupportRecip);

                // BRAKING-GATED RESPAN (v0.14). While braking, an
                // upward release retarget re-derives the span entry-
                // style (relTarget-prev) instead of incrementally. The
                // brake's own retargets — a settle handing its floor
                // back to lookaheadX — are exactly the large upward
                // jumps where the grown incremental span and the
                // unchanged remaining distance go inconsistent and the
                // accelerating anchor would lunge past the target by
                // the already-traversed distance. Re-deriving the span
                // restores consistency by construction, so the reach
                // rule (Step 5) can pick freely: measured, this nearly
                // halves the pre-respan residual flips while returning
                // 0.016-0.024 of mean envelope. Outside braking the
                // incremental span stands — per-jitter respans on
                // plain noise measurably cost level — and with the
                // brake disabled the gate strips the respan entirely.
                totalStep = select2(releasing,
                    min(select2(prevTotalStep>=0, incrementalTotalStep, entryTotalStep),
                        supportTotalStep),
                    max(select2((prevTotalStep<=0)|((relTarget>relTarget')&braking), incrementalTotalStep, entryTotalStep),
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
                // The last test moved off prevPhase == 0: the rise leg of a
                // turnaround ends at phase 0.0 exactly, and the descent
                // re-enters through it.
                sameTarget = (relTarget==relTarget');
                sameDirection = (releasing==releasing')&(attacking==attacking');
                undisturbed = (prev==prevExpected);
                trusted = sameTarget&sameDirection&undisturbed&(prevTotalStep!=0);

                // --- Step 5: velocity matching (for untrusted samples) ---
                //
                // Since the emitted delta is an exact curve increment, the
                // measured speed is the MEAN derivative over the last step.
                // The inverse-derivative functions invert the POINT
                // derivative, so they return (approximately) the midpoint of
                // that step. Adding halfStep converts the recovered midpoint
                // into the current position — on either side of the apex:
                // the rise phase advances by +step like any other, so the
                // recovered |midpoint| is current depth + halfStep and the
                // correction is additive after negation too.
                //
                // The correction only makes sense when a previous interval
                // in the same direction actually exists; the gate is the
                // same in-transition sentinel as Step 4.
                matchCorr = halfStep*(sameDirection&(prevTotalStep!=0));
                speedRatio = prevSpeed/(totalStep*invCurveScale*step+(1-active)*1e-30);

                // Shared inverse-derivative inversion (v0.7): one sqrt
                // and one rational denominator serve every anchor below.
                matchD = max(0, min(speedRatio, maxDerivVal));
                partD = inverseDerivativePart(shape, matchD);
                denD = 2*(shape*matchD+1);
                topReleaseD = (1+partD)/denD;
                bottomReleaseD = (1-partD)/denD;

                // gonnaMakeIt projects the landing from the decelerating
                // anchor and picks by whether the trajectory still reaches
                // the target. The preference is deliberate: an anchor that
                // reaches can at worst overshoot, which the clamps and the
                // landing guard catch; an anchor that strands short
                // arrives LATE, and on the attack side late is the one
                // thing the lookahead contract cannot absorb.
                //
                // The projection needs the curve's fraction-done at the
                // late anchor. Two ways to get it, one compile-time switch
                // (exactBranchPick):
                //
                //   default: evaluate the per-direction tail cubics fitted
                //   at slider rate (see the tail fit above). The anchor
                //   always lies on [peakPhase, 1], so the fit covers
                //   everything the projection can ask for. The pick can
                //   differ from the exact one only when the projected
                //   landing sits within the fit error of the target
                //   (2e-4 of the span for the attack, 1.3e-2 for the
                //   release) — where both picks land within a few guard
                //   sizes and either is absorbed: a snap or one small
                //   fresh transition. The attack's duck deadline is
                //   protected to 2e-4 of the span; the release has no
                //   deadline.
                //
                //   exact (v0.4): evaluate the curve itself — the
                //   log+atan pair, the single most expensive thing in the
                //   loop — every sample.
                lateAnchor = select2(releasing,
                    1-bottomReleaseD,
                    topReleaseD)+matchCorr:min(1);

                gonnaDo(phase) = (1-fdOf(phase))*totalStep;
                projected = gonnaDo(lateAnchor)+prev;

                uLate = (lateAnchor-select2(releasing, attackPeakPhase, releasePeakPhase))
                    *select2(releasing, 1/(1-attackPeakPhase), 1/(1-releasePeakPhase))
                    :max(0):min(1);
                tailPiece = uLate>tailSplit;
                tailT = (uLate-tailPiece*tailSplit)
                    *select2(tailPiece, 1/tailSplit, 1/(1-tailSplit));
                tcSel(tc) = select2(releasing,
                    select2(tailPiece,
                        tc(attLoF(0), attLoF(1), attLoF(2), attLoF(3)),
                        tc(attHiF(0), attHiF(1), attHiF(2), attHiF(3))),
                    select2(tailPiece,
                        tc(relLoF(0), relLoF(1), relLoF(2), relLoF(3)),
                        tc(relHiF(0), relHiF(1), relHiF(2), relHiF(3))));
                fdLate = tcSel(tc0)+tailT*(tcSel(tc1)+tailT*(tcSel(tc2)+tailT*tcSel(tc3)));
                projectedFast = (1-fdLate)*totalStep+prev;

                gonnaMakeIt = select2(exactBranchPick,
                    projectedFast>lookaheadX,
                    projected>lookaheadX);

                // THE REACH RULE, RESTORED (v0.14). The release picks
                // the accelerating anchor unless the decelerating
                // anchor still reaches the target — v0.4's rule, same
                // as the attack's. v0.9 had replaced it with a strand-
                // versus-overshoot comparison to duck a real lunge:
                // after an upward retarget the incremental span has
                // grown but the remaining distance has not, so the
                // accelerating anchor overshoots by the already-
                // traversed distance. That patch treated the symptom
                // and billed every noisy rise for it — its overshoot
                // estimate overstates systematically late in flight,
                // so dense upward retargets kept re-anchoring releases
                // onto the decelerating tail, and the envelope rounded
                // off and stalled far below its target: 0.03-0.04 of
                // mean envelope across the noisy battery, the
                // premature-deceleration stall on the scope. The
                // lunge's root — the span/position inconsistency — is
                // now fixed at the source while braking (the respan,
                // Step 3); outside braking the surviving lunges
                // measure at v0.4 parity (max ~0.004 of full scale,
                // zero deep corners): the cost the reference behavior
                // always paid, caught by the ordinary attack-back.
                forwardPos = select2(releasing,
                    // Attack branch
                    select2(gonnaMakeIt,
                        1-bottomReleaseD,
                        1-topReleaseD),
                    // Release branch
                    select2(gonnaMakeIt,
                        bottomReleaseD,
                        topReleaseD))+matchCorr:max(0):min(1-step);

                matchedPos = forwardPos;

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
                // quadrature on the cheap rational derivative (see v0.4
                // note B; the derivative of fdOf is the same expression for
                // both directions — the attack flip's two sign changes
                // cancel in the chain rule). The select keeps attack
                // increments exact; release increments take the larger
                // of the brake's two allowances (see Step 1).
                glA = 0.21132486540518713;
                glB = 0.78867513459481287;
                derivOf(p) = derivativeBaseRelease(shape,
                    select2(releasing, 1-p, p))*invCurveScale;
                // RISE WHILE YOU CAN (v0.14). Applies to the ABOVE-
                // floor branch only — a release already over the
                // coming minimum, the one that must reach zero speed
                // by the bind. It is exempt from the fade when its
                // remaining flight fits inside a TRUSTED countdown: it
                // lands — zero speed — before the bind, so the bump
                // toward the present target is taken in full and the
                // attack departs from a standstill at the target
                // itself. One that cannot land in time is faded
                // linearly to zero exactly at the bind: a partial
                // bump, still kink-free. Below the floor nothing
                // changes (settle onto steady floors, ride slopes;
                // such a release never needs zero bind speed), and
                // only countdown-invalid standstill births remain
                // held. Fading below-floor releases too — tried and
                // reverted — reimports the v0.13 level loss on dense
                // material: perpetual gaps keep the ramp saturated and
                // the fade freezes releases that threaten nothing.
                willLand = (1-anchor)*durSamples
                    <= (1-brakeRamp)*brake_samples;
                fading = releasing&braking&(prev>=attackMinAhead)
                    &(1-(gapFresh&willLand));
                delta = totalStep*step*0.5
                    *(derivOf(anchor+glA*step)+derivOf(anchor+glB*step))
                    *(1-brakeRamp*fading);

                // The value this loop wants to write. Stored as state so the
                // next sample can detect outside modification (see Step 4).
                expected = prev+delta;

                // Completion: a transition ends where its phase runs
                // out; the halfStep slack absorbs accumulation rounding
                // in the final phase value (the phase advances in whole
                // steps, so the test is unambiguous).
                //
                // Landing: snap to the target — but only across a
                // residual the snap was designed for (quadrature epsilon,
                // float epsilon). guardSize is one peak-speed sample.
                // Under mid-transition retargets the phase can complete
                // with a real gap left (see v0.4 note A); snapping across
                // that is a teleport, so instead the state resets WITHOUT
                // snapping and the next sample starts a fresh shaped
                // transition over the leftover distance, entered through
                // velocity matching so it continues from the current
                // speed.
                //
                // The output never goes above the raw (delayed) input and
                // never below the attack target. The lower clamp is the
                // required catcher for supportTotalStep's deliberate
                // overshoot (Step 3); it is inert outside of attack, since
                // min(prev, lookaheadX) equals prev while releasing or
                // idle. It cannot conflict
                // with the upper clamp because lookaheadX <= delayedX by
                // construction (the sliding-min window contains the
                // delayed sample).
                phaseDone = active&(newPhaseRaw>=(1-halfStep));
                guardSize = abs(totalStep)*maxDerivVal*invCurveScale*step;
                landed = phaseDone&(abs(relTarget-expected)<=guardSize);

                result = select2(landed,
                    min(expected, delayedX),
                    min(relTarget, delayedX)):max(min(prev, lookaheadX));

                // With the state cleared: constant target -> idle as
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
//    Progress through the current transition, on [0, 1]. Advances by
//    `step` per sample; trusted between discontinuities, re-derived
//    from the measured speed (velocity matching) across them.
//
//  braking (v0.12, re-clocked v0.13, finished v0.14):
//    The future attack target (attackMinAhead = lookaheadX read
//    brake_samples ahead) is below the present one; a bind is exactly
//    brake_samples away when this turns true. A release above the
//    coming floor is faded out by 1-brakeRamp, a linear ramp that
//    saturates exactly at the bind, so its attack enters from a
//    standstill; a release below a STEADY coming floor (one deepest
//    sample reigning: attackMinAhead bit-constant) retargets to the
//    floor itself and runs an ordinary shaped release into it,
//    landing and parking (v0.14); drifting minima — slopes —
//    are ridden at the attack window. ABOVE the floor (v0.14), a
//    release whose remaining flight fits a trusted countdown runs
//    free and lands before the bind; one that cannot is faded to
//    zero exactly at it. Standstill births under a pending
//    at-or-below duck are held (freshHold): recoveries shorter than
//    the brake window are suppressed, release-hold style.
//
//  turnaround (v0.5–v0.11, retired):
//    The mirrored-curve reversal scheduled by an arrival countdown.
//    Replaced by the brake; see the v0.12 paragraph at the top.
//
//  arrival (v0.5–v0.11, retired):
//    The countdown that scheduled the turnaround's just-in-time flip.
//    Gone with it: attacks start the moment the attack window binds.
//
//  holdCeil (v0.10–v0.11, retired):
//    The latched pre-episode minimum that capped the rise apex. The
//    brake never lifts the envelope above anything in flight.
//
//  held (v0.5–v0.11, retired):
//    The just-in-time entry gate that delayed wanted attacks while
//    schedule budget remained. Attacks are no longer gated.
//
//  totalStep:
//    The full amplitude span (in linear gain) of the current transition
//    leg. Maintained incrementally: held exactly while the target is
//    static, moved by exactly (lookaheadX - lookaheadX') when it
//    isn't, recomputed on entry as (lookaheadX - prev), and never
//    allowed to shrink past the span whose curve can still carry the
//    current speed (supportTotalStep, two-sided on the attack).
//
//  step:
//    Phase increment per sample = 1 / (duration_in_samples).
//
//  fdOf(p):
//    Fraction of the transition completed at phase p, for both directions.
//    Position on the curve = transitionStart + totalStep*fdOf(p). The
//    per-sample delta integrates fdOf' over one step with 2-point
//    Gauss-Legendre quadrature; the ~1e-10 per-transition residual is
//    absorbed by the landing snap.
//
//  velocity matching:
//    Across a discontinuity (retarget, direction flip, outside
//    modification), find the phase on the (possibly new) curve whose
//    derivative equals the envelope's current speed, then continue from
//    there. A speed whose sign opposes the span clamps to the
//    zero-derivative anchor — the leg's start; e.g. an attack entered
//    from a braked, still-creeping release. With exact increments the
//    measured speed is a mean over the last step, hence the halfStep
//    matchCorr.
//
//  gonnaMakeIt:
//    "If we re-anchor on the decelerating branch now, will the resulting
//    trajectory still reach the target?" Picks the accelerating vs
//    decelerating branch of the inverse derivative, preferring the one
//    that reaches (overshoot is clampable; stranding short is late). By
//    default the projection evaluates slider-rate-fitted tail cubics;
//    exactBranchPick = 1 evaluates the true curve (the log+atan pair)
//    every sample instead.
//
//  landing:
//    When a leg's phase completes, the transition state is fully reset.
//    A forward leg snaps to the target only if the residual is at most
//    one peak-speed sample (quadrature/float epsilon territory); larger
//    residuals are NOT snapped (teleport) — the leftover gets a fresh
//    shaped transition that velocity-matches into the current speed. The
//    phase-completion test carries halfStep of slack because the final
//    phase value is produced by accumulation and can round off-grid. An
//    attack starts the sample the attack window binds and takes
//    att_samples+1 — landing exactly when delayedX presents the
//    transient.
//
//  lookaheadX:
//    The sliding minimum of the input over the next att_samples (the
//    attack window): the envelope's target. attackMinAhead is the
//    same minimum one brake window earlier in its pipeline — the
//    future target — and the brake clock compares the two.
//
//  AUC compensation:
//    (Omitted from this core file for clarity.) Adjusts the attack/release
//    duration so that changing the shape doesn't change the average level
//    of the envelope — sharper shapes get proportionally shorter durations.
