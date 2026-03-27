declare name "shapedSmoother";
declare version "0.2";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

process = test2@att_samples, shapedSmoother(test2);

shapedSmoother(x) = (lookaheadX:env~(_, _, _))//:(_, _, !, _)// :(_, !, !)
    with {
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR);
        //@half_att_samples;
        env(prev, prevIndex, prevTotalStep, lookaheadX) = result, newIndex, totalStep
            with {
                shape = shapeMap(select2(releasing,
                    hslider("attack shape", 0, 0, 1, 0.001),
                    hslider("release shape", 0, 0, 1, 0.001)));
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active = attacking|releasing;

                // Time base depends on direction
                activeTime = select2(releasing, att, rel);
                totalNRSteps = (activeTime*ma.SR):max(1);
                step = 1/totalNRSteps;

                // Total step: min for attack (negative), max for release (positive)

                todo = lookaheadX-prev;

                remainingCurve = 1-select2(releasing',
                    cheapCurveAttack(shape, prevIndex/(totalNRSteps')),
                    cheapCurveRelease(shape, prevIndex/(totalNRSteps')));
                // totalStep = todo/max(1/totalNRSteps, remainingCurve);
                totalStep = todo/remainingCurve;

                old_totalStep = select2(releasing,
                    (lookaheadX-prev):min(prevTotalStep),
                    (lookaheadX-prev):max(prevTotalStep))*active;

                // Convert integer index to phase for curve calculations
                prevPhase = prevIndex/(totalNRSteps');

                // Attack curve: slow start, fast end. Release curve: fast start, slow end.
                prevSpeed = prevTotalStep*select2(releasing,
                    derivativeAttack(shape, prevPhase),
                    derivativeRelease(shape, prevPhase));
                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0, min(speedRatio, maxDerivativeAttack(shape)));

                // Convert new index to phase for speed calculation
                newPhase = newIndex/totalNRSteps;
                // newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                speed = totalStep*select2(releasing,
                    derivativeAttack(shape, newPhase),
                    derivativeRelease(shape, newPhase));

                gonnaDo(phase) = (1-select2(releasing,
                    cheapCurveAttack(shape, phase),
                    cheapCurveRelease(shape, phase)))*totalStep;

                // TODO: find proper fix for wrong values of gonnaMakeIt, these "+step*0.x" are just a workaround
                projected = gonnaDo(select2(releasing,
                    inverseDerivativeBottomAttack(shape, clampedRatio)+(step*0.3),
                    inverseDerivativeTopRelease(shape, clampedRatio)+(step*0.5))// +(step*hslider("step mult", 0, 0, 1, 0.1))
                )+prev;
                gonnaMakeIt = (projected>lookaheadX);

                phaseAtMatchingSpeed = select2(releasing,
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomAttack(shape, clampedRatio),
                        inverseDerivativeTopAttack(shape, clampedRatio)),
                    select2(gonnaMakeIt,
                        inverseDerivativeBottomRelease(shape, clampedRatio),
                        inverseDerivativeTopRelease(shape, clampedRatio)));

                newIndex = (phaseAtMatchingSpeed*totalNRSteps+1):min((totalNRSteps-1)):max(1)*active;

                shaped = (speed*step)+prev;
                result = min(shaped, x@att_samples);
            };
    };
// Parameters
maxHold = 0.05;
// maxSR = 192000;
maxSR = 48000;
att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.001)/1000;
att_samples = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
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
derivativeBaseRelease(c, x) = x*(1-x)/(c*pow(x, 2)+1-c);
derivativeBaseAttack(c, x) = x*(1-x)/(c*pow(x, 2)+1-(2*c*x));
derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x) = derivativeBaseAttack(c, x)/curveScale(c);
peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
maxDerivativeAttack(c) = derivativeAttack(c, peakPhaseAttack(c));
inverseDerivativePart(c, x) = sqrt(max(0, 1-4*(curveScale(c)*x-c*curveScale(c)*x)*(c*curveScale(c)*x+1)));
inverseDerivativeTopRelease(c, x) = (1+inverseDerivativePart(c, x))/(2*(c*curveScale(c)*x+1));
inverseDerivativeBottomRelease(c, x) = (1-inverseDerivativePart(c, x))/(2*(c*curveScale(c)*x+1));
inverseDerivativeTopAttack(c, x) = inverseDerivativeTopRelease(c, x)*-1+1;
inverseDerivativeBottomAttack(c, x) = inverseDerivativeBottomRelease(c, x)*-1+1;
shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));
