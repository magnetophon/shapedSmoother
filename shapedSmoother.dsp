declare name "shapedSmoother";
declare version "0.3.2";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

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

process = test2@(att_samples+brake_samples), shapedSmoother(test2);
shapedSmoother(x) = lookaheadX:env~(_, _, _)
    with {
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHoldSamples)@brake_samples;

        // brake: covers an additional brake_samples ahead, used purely for
        // early attack detection.
        lookahead_brake = x:ba.slidingMin(att_samples+brake_samples+1, 1+(maxHold+maxBrake)*maxSR);
        // ---- Hoisted per-shape constants (slider-rate) ----
        // These are pure functions of the shape sliders. Faust's auto-lift
        // doesn't push through `select2(releasing, sliderA, sliderB)` so in
        // v0.3 every `cheapCurveBase`/`curveScale`/`maxDerivativeBaseAttack`
        // call on `shape` fires audio-rate. Precompute them per-shape.
        attackShape = shapeMap(hslider("attack shape", 0, 0, 1, 0.001));
        releaseShape = shapeMap(hslider("release shape", 0, 0, 1, 0.001));

        attackInvCurveScale = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackZero = cheapCurveBase(attackShape, 0);
        releaseZero = cheapCurveBase(releaseShape, 0);

        attackMaxDerivBase = maxDerivativeBaseAttack(attackShape);
        releaseMaxDerivBase = maxDerivativeBaseAttack(releaseShape);

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

                prevSpeed = prevTotalStep*select2(releasing,
                    derivativeBaseAttack(shape, prevPhase),
                    derivativeBaseRelease(shape, prevPhase));
                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0,
                    min(speedRatio,
                        maxDerivVal));

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                gonnaDo(phase) = (1-select2(releasing,
                    cheapCurveAttackS(shape, invCurveScale, zeroVal, phase),
                    cheapCurveReleaseS(shape, invCurveScale, zeroVal, phase)))*totalStep;

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio)))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                phaseAtMatchingSpeed = select2(releasing,
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        inverseDerivativeTopRelease(shape, clampedRatio)));

                speed = totalStep*select2(releasing,
                    derivativeAttackS(shape, invCurveScale, newPhase),
                    derivativeReleaseS(shape, invCurveScale, newPhase));

                braking = (prev>lookahead_brake)&(prev<=lookaheadX);
                brakeRamp = ((_+(1/brake_samples))*braking)~_:min(1);
                brakeMult = 1-brakeRamp;
                delta = speed*step*brakeMult;

                result = min(prev+delta, x@(att_samples+brake_samples));
            };
    };

// Parameters
maxHold = 0.05;
maxBrake = 0.05;
maxSR = 48000;

maxHoldSamples = maxHold*maxSR;
maxBrakeSamples = maxBrake*maxSR;
maxTotalSamples = (maxHold+maxBrake)*maxSR;

att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01)/1000;
att_samples = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
rel = hslider("rel[scale:log]", 0.05*1000, 1, 5000, 0.1)/1000;
rel_samples = rel*ma.SR:max(1);

// Brake window length per spec:
//   = 50 ms by default,
//   = max(att, rel) if both att AND rel are shorter than 50 ms.
// (min(50ms, max(att,rel)) gives exactly this.)
brake_samples = min(maxBrakeSamples, max(att_samples, rel_samples)):max(1);

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
