declare name "shapedSmoother";
declare version "0.3";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

// lookahead array 
// LA1diff*0.5 
// LA2
// 
// check if:
// - we need to go down steeper than we are going now -> 
linearProj = prev+(n*prevDelta);
steeperDown = linearProj>lookaheadBrake;
keepSpeed = linearProj==lookaheadBrake;
lessSteepDown = linearProj<lookaheadBrake;
// keep count of how many samples the LA has been the same
// keep count of how many samples the LA has gone down
// 
// 
//
// dirToPoint(dir, point, nr) =

process = test2@att_samples, att_env(test2);
att_env(x) = env~(_, _, _)//
//:(_, _, !, _)
:(_, _, (_*att_samples), !)
    with {
        env(prevSmoother, prevDelta, prevTargetRate) = smoother, delta, targetRate, targetSameCounter//, LA
            with {
                remainingSteps = max(1, (1-targetSameCounter)*totalNRSteps);
                targetRate = ((LA-prevSmoother)/remainingSteps);
                delta = (prevDelta, targetRate):it.interpolate_linear(blendFactor);
                blendPow = hslider("blend curve", 3, 1, 8, 0.1);
                // blendFactor = pow(targetSameCounter, blendPow);
                // blendFactor = cheapCurveAttack(shape, targetSameCounter);
                blendK = hslider("blend curve", 5, 1, 10, 0.1);
                blendFactor = ((exp(blendK*targetSameCounter)-1)/max(1e-10, exp(blendK)-1));

                // delta = ((derivativeAttack(shape, targetSameCounter)), (LA-prevSmoother)):it.interpolate_linear(targetSameCounter);

                // delta = (prevDelta, targetRate*targetSameCounter):it.interpolate_linear(targetSameCounter);
                // delta = (prevDelta, targetRate):it.interpolate_linear(targetSameCounter:cheapCurveAttack(shape));
                smoother = (prevSmoother+delta):min(x@att_samples);
                targetSameCounter = ((_+step)*same:min(1))~_;
                LA = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR);
                // delta = // delta = ((derivativeAttack(shape, targetSameCounter)));
                // ((prevSmoother'-prevSmoother), (LA-prevSmoother)):it.interpolate_linear(targetSameCounter);
                totalNRSteps = att_samples;
                step = 1/totalNRSteps;
                same = LA==LA';
                // lookaheadX==lookaheadBrake;
                lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR)@half_att_samples;
                lookaheadBrake = x:ba.slidingMin(brake_samples+1, 1+1.5*maxHold*maxSR);

                shape = shapeMap(hslider("attack shape", 0, 0, 1, 0.001));
                active = prev!=LA;
            };
    };

// lookaheadX(x) = x:ba.slidingMin(32+1, 1+32)@half_att_samples;
// lookaheadBrake(x) = x:ba.slidingMin(hslider("LA", 0, 0, 32, 1)+1, 32);

// test2@brake_samples, shapedSmoother(test2);

shapedSmoother(x) = (x:env~(_, _, _))//
//:(_, !, !, _, _)
    with {
        env(prev, prevPhase, prevTotalStep, x) = result, newPhase, totalStep//, delta//, lookaheadX, lookaheadBrake
            with {
                lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR)@half_att_samples;
                lookaheadBrake = x:ba.slidingMin(brake_samples+1, 1+1.5*maxHold*maxSR);
                shape = shapeMap(select2(releasing,
                    hslider("attack shape", 0, 0, 1, 0.001),
                    hslider("release shape", 0, 0, 1, 0.001)));
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                activeTime = select2(releasing, att, rel);
                totalNRSteps = (activeTime*ma.SR):max(1);
                step = 1/totalNRSteps;

                // make delta dependent on how much GR we have to do until the full lookahead, so it self corrects

                todo = lookaheadX-prev;

                remainingCurve = 1-select2(releasing,
                    cheapCurveAttack(shape, prevPhase),
                    cheapCurveRelease(shape, prevPhase));
                new_totalStep = todo/max(1/totalNRSteps, remainingCurve);
                // totalStep = todo/remainingCurve;

                // TODO: clamp the totalStep so that the speed at peakPhaseAttack matches the delta 
                // totalStep = (todo/remainingCurve)<:select2(releasing,
                totalStep = (todo)<:select2(releasing,
                    min(prevTotalStep),
                    max(select2(lookaheadBrake<lookaheadX, prevTotalStep, 0-ma.INFINITY)));
                prevSpeed = prevTotalStep*select2(releasing,
                    derivativeBaseAttack(shape, prevPhase),
                    derivativeBaseRelease(shape, prevPhase));
                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0,
                    min(speedRatio,
                        maxDerivativeBaseAttack(shape)));

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                gonnaDo(phase) = (1-select2(releasing,
                    cheapCurveAttack(shape, phase),
                    cheapCurveRelease(shape, phase)))*totalStep;

                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio)))+prev;
                projected2 = gonnaDo(select2(releasing,
                    inverseDerivativeTopAttack(shape, clampedRatio),
                    inverseDerivativeBottomRelease(shape, clampedRatio)))+prev;
                gonnaMakeIt = (projected>lookaheadX);

                phaseAtMatchingSpeed = select2(releasing,
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        inverseDerivativeTopRelease(shape, clampedRatio)));

                // Final output still uses the scaled derivative (needed for
                // the actual envelope value computation)
                speed = totalStep*select2(releasing,
                    derivativeAttack(shape, newPhase),
                    derivativeRelease(shape, newPhase));

                delta = speed*step;
                result = min(prev+delta, x@brake_samples);
            };
    };

// Parameters
maxHold = 0.05;
maxSR = 48000;
att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.01)/1000;
half_att_samples = (0.5*att*ma.SR):max(1);
att_samples = 2*half_att_samples;
brake_samples = 3*half_att_samples;
rel = hslider("rel[scale:log]", 0.05*1000, 1, 5000, 0.1)/1000;

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

// Base (unscaled) derivatives — these are the raw x(1-x)/denominator without
// dividing by curveScale.
derivativeBaseRelease(c, x) = x*(1-x)/(c*pow(x, 2)+1-c);
derivativeBaseAttack(c, x) = x*(1-x)/(c*pow(x, 2)+1-(2*c*x));

// Scaled derivatives (divided by curveScale) — used only where the actual
// envelope value is needed, not in the inverse path.
derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);

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
