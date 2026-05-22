declare name "shapedSmoother";
declare version "0.4";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

process = test2@att_samples, shapedSmoother(test2);

shapedSmoother(x) = (lookaheadX:env~(_, _, _))
    with {
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR);

        // ---- Control-rate per-shape constants -----------------------------
        // The original code routes `shape` through select2(releasing, ...),
        // which makes it a sample-rate signal (releasing depends on the
        // audio-rate feedback `prev`). Every shapeMap / curveScale / peak /
        // maxDeriv call downstream of that select is therefore forced to run
        // at audio rate.
        //
        // By precomputing the heavy math per-shape (control rate) and
        // selecting only the scalar result audio-rate, exp / pow / log /
        // atan / sqrt all move to the slider-change path.
        attackShape  = shapeMap(hslider("attack shape",  0, 0, 1, 0.001));
        releaseShape = shapeMap(hslider("release shape", 0, 0, 1, 0.001));

        // 1/curveScale instead of curveScale: turns the per-sample division
        // into a per-sample multiply.
        attackInvCurveScale  = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        // cheapCurveBase(c, 0) shows up inside every cheapCurveRelease call.
        // Hoist it.
        attackZero  = cheapCurveBase(attackShape,  0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        // Unscaled peak derivative — used for clamping speedRatio. Due to the
        // horizontal-flip symmetry the peak value of the release curve
        // equals the peak value of the attack curve, but the peak is evaluated
        // at the current per-mode shape value either way.
        attackMaxDeriv  = maxDerivativeBaseAttack(attackShape);
        releaseMaxDeriv = maxDerivativeBaseAttack(releaseShape);

        // Pre-split 1/(activeTime*SR) per mode for the same reason.
        attNRSteps = (att*ma.SR):max(1);
        relNRSteps = (rel*ma.SR):max(1);
        attStep    = 1/attNRSteps;
        relStep    = 1/relNRSteps;

        env(prev, prevPhase, prevTotalStep, lookaheadX) = result, newPhase, totalStep
            with {
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active    = attacking|releasing;

                // Audio-rate scalars: pick the right per-shape constant. The
                // heavy math already happened above.
                shape         = select2(releasing, attackShape,         releaseShape);
                invCurveScale = select2(releasing, attackInvCurveScale, releaseInvCurveScale);
                zeroVal       = select2(releasing, attackZero,          releaseZero);
                maxDerivVal   = select2(releasing, attackMaxDeriv,      releaseMaxDeriv);
                step          = select2(releasing, attStep,             relStep);

                todo = lookaheadX-prev;
                totalStep = select2(releasing,
                    todo:min(prevTotalStep),
                    todo:max(prevTotalStep)) * active;

                // ---- Previous-sample derivative -------------------------
                // derivativeBaseAttack(c, x) == derivativeBase(c, 1-x), so
                // both branches collapse into one phase-select + one call.
                prevY     = select2(releasing, 1-prevPhase, prevPhase);
                prevSpeed = prevTotalStep * derivativeBase(shape, prevY);

                speedRatio   = prevSpeed/totalStep;
                clampedRatio = max(0, min(speedRatio, maxDerivVal));

                // ---- Inverse derivative ---------------------------------
                // All four originals (top/bottom × attack/release) share
                //   inverseDerivativePart(c, D)  = sqrt(...)
                //   1/(2*(c*D+1))
                // Compute these shared pieces once. Attack variants are
                // simply (1 - their release counterpart) and are folded into
                // the consumers below rather than computed eagerly.
                cD1     = shape*clampedRatio + 1;
                invPart = sqrt(max(0, 1 - 4*clampedRatio*(1-shape)*cD1));
                invDen2 = 1/(2*cD1);
                topRelease    = (1 + invPart)*invDen2;
                bottomRelease = (1 - invPart)*invDen2;

                // ---- Projected target ("gonnaDo + prev") ----------------
                // Working through cheapCurveAttack(c, 1-bottomRelease):
                //   1 - cheapCurveAttack(c, 1-b) == cheapCurveRelease(c, b)
                // so the attack and release branches both reduce to one
                // cheapCurveRelease call at a selected phase, with a flip on
                // release. That eliminates one cheapCurveBase evaluation per
                // sample (a log + an atan).
                gonnaDoPhase  = select2(releasing, bottomRelease, topRelease);
                gonnaDoCurve  = (cheapCurveBase(shape, gonnaDoPhase) - zeroVal) * invCurveScale;
                gonnaDoFactor = select2(releasing, gonnaDoCurve, 1 - gonnaDoCurve);
                projected     = gonnaDoFactor*totalStep + prev;
                gonnaMakeIt   = (projected>lookaheadX);

                // ---- Phase at matching speed ----------------------------
                // The 2x2 (releasing × gonnaMakeIt) reduces to one inner
                // select + an outer flip, because attack values are just
                // (1 - their release counterparts).
                inner = select2(gonnaMakeIt, bottomRelease, topRelease);
                phaseAtMatchingSpeed = select2(releasing, 1 - inner, inner);

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                // ---- Final speed (scaled derivative at newPhase) --------
                newY  = select2(releasing, 1-newPhase, newPhase);
                speed = totalStep * derivativeBase(shape, newY) * invCurveScale;

                result = min(prev+speed*step, x@att_samples);
            };
    };

// Parameters
maxHold = 0.05;
maxSR   = 48000;
att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01)/1000;
att_samples      = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
rel = hslider("rel[scale:log]", 0.05*1000, 1, 5000, 0.1)/1000;

// Test signals
test2 = it.interpolate_linear(hslider("noise level", 0, 0, 1, 0.001),
    (loop~_),
    no.lfnoise(hslider("noise rate", 42, 1, 1000, 1)))
    with {
        loop(prev) = no.lfnoise0(blockscale*(abs(prev*69)%9:pow(0.75)*5+1));
        blockscale = hslider("blockscale", 1, 0.01, 10, 0.01);
    };
testSig = os.lf_sawpos(0.5)<:(((_>0.25)*hslider("step1", 0.75, -1, 1, 0.001))+((_>0.5)*hslider("step2", 0.125, -1, 1, 0.001)));

// *************************************** the NEW curves ******************************
// nuchi's curves — same math as before, reorganised so attack and release
// share work (they're horizontal flips of each other):
//   https://www.wolframalpha.com/input?i=antiderivative+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29
//   https://www.wolframalpha.com/input?i=inverse+function+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29+%2F+C
//   https://www.desmos.com/calculator/9wtfhymvr0
//   scaled c:        https://www.desmos.com/calculator/dynyjjkuli
//   horizontal flip: https://www.desmos.com/calculator/bsr8cdn21v

cheapCurveBase(c, x) = (log(c*(x*x-1)+1)+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);
curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x)  = 1 - cheapCurveRelease(c, 1-x);

// Unified unscaled derivative: derivativeBase(c, y) with y=x for release and
// y=1-x for attack. Algebraic proof:
//   derivativeBaseAttack(c, x) = x(1-x) / (c*x^2 + 1 - 2*c*x)
//                              = y(1-y) / (c*y^2 + 1 - c)   with y = 1-x
//                              = derivativeBase(c, 1-x)
// This lets the env() collapse two function calls into one.
derivativeBase(c, y) = y*(1-y)/(c*y*y + 1 - c);

derivativeBaseRelease(c, x) = derivativeBase(c, x);
derivativeBaseAttack(c, x)  = derivativeBase(c, 1-x);

derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x)  = derivativeBaseAttack(c, x)/curveScale(c);

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));
maxDerivativeAttack(c)     = derivativeAttack(c, peakPhaseAttack(c));

// Inverse-derivative pieces kept for reference; the hot path in env() inlines
// the shared sqrt + 1/(2(cD+1)) directly to avoid recomputing them across
// four near-identical functions.
inverseDerivativePart(c, D)          = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));
inverseDerivativeTopRelease(c, D)    = (1+inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeBottomRelease(c, D) = (1-inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeTopAttack(c, D)     = 1-inverseDerivativeTopRelease(c, D);
inverseDerivativeBottomAttack(c, D)  = 1-inverseDerivativeBottomRelease(c, D);

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));
