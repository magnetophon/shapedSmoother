declare name "shapedSmoother";
declare version "0.3";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

process = test2@att_samples, shapedSmoother(test2);

shapedSmoother(x) = x//
:releaseEnv~(_, _, _)//
//:(_, !, !)//
:(_, _, !)//
:attackEnv~(_, _, _):(_, _, !, _)
    with {
        // Release envelope: smooths upward movement, no lookahead.
        // When not releasing, tracks input instantly (passes through drops).
        releaseEnv(prev, prevPhase, prevTotalStep, x) = result, newPhase, totalStep//
            with {
                shape = shapeMap(hslider("release shape", 0, 0, 1, 0.001));
                releasing = (prev<x)&(rel>0);
                active = releasing;

                totalNRSteps = (rel*ma.SR):max(1);
                step = 1/totalNRSteps;

                totalStep = ((x-prev):max(prevTotalStep))*active;

                prevSpeed = prevTotalStep*derivativeRelease(shape, prevPhase);
                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0, min(speedRatio, maxDerivativeAttack(shape)));

                speed = totalStep*derivativeRelease(shape, newPhase);

                gonnaDo(phase) = (1-cheapCurveRelease(shape, phase))*totalStep;

                projected = gonnaDo(inverseDerivativeTopRelease(shape, clampedRatio)+(step*0.5))+prev;
                gonnaMakeIt = (projected>x);

                phaseAtMatchingSpeed = select2(gonnaMakeIt,
                    inverseDerivativeBottomRelease(shape, clampedRatio),
                    inverseDerivativeTopRelease(shape, clampedRatio));

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

                shaped = (speed*step)+prev;
                // When not active, pass through x so prev tracks input (instant drop-following).
                // min(_, x) clamps release overshoot.
                result = min(select2(active, x, shaped), x);
            };

        // Attack envelope: smooths downward movement, with lookahead.
        // When not attacking, tracks input instantly (passes through rises).
        attackEnv(prev, prevPhase, prevTotalStep, y, releasePhase) = result, newPhase, totalStep//
        , releasePhase@att_samples*(active==0)//
        // , result<=lookaheadY
            with {
                lookaheadY = y:ba.slidingMin(att_samples+1, 1+maxHold*maxSR);

                shape = shapeMap(hslider("attack shape", 0, 0, 1, 0.001));
                attacking = (prev>=lookaheadY)&(att>0);
                active = attacking;
                holding = (attacking==0)&((y@att_samples-step)>prev)&(prevPhase>0);

                totalNRSteps = (att*ma.SR):max(1);
                step = 1/totalNRSteps;

                totalStep = ((lookaheadY-prev):min(prevTotalStep))*active;

                prevSpeed = prevTotalStep*derivativeAttack(shape, prevPhase);
                speedRatio = prevSpeed/totalStep;
                clampedRatio = max(0, min(speedRatio, maxDerivativeAttack(shape)));

                speed = totalStep*derivativeAttack(shape, newPhase);

                gonnaDo(phase) = (1-cheapCurveAttack(shape, phase))*totalStep;

                projected = gonnaDo(inverseDerivativeBottomAttack(shape, clampedRatio)+(step*0.3))+prev;
                gonnaMakeIt = (projected>lookaheadY);

                phaseAtMatchingSpeed = select2(gonnaMakeIt,
                    inverseDerivativeBottomAttack(shape, clampedRatio),
                    inverseDerivativeTopAttack(shape, clampedRatio));

                // newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;
                newPhase = select2(holding,
                    (phaseAtMatchingSpeed+step):min(1-step):max(step)*active,
                    prevPhase);

                shaped = (speed*step)+prev;
                // When not active, pass through y@att_samples so prev tracks input (instant rise-following).
                // min(_, y@att_samples) clamps attack overshoot.
                // result = min(select2(active, y@att_samples, shaped), y@att_samples);
                // result = select2(active, y@att_samples, shaped);
                result = select2(active, select2(holding, y@att_samples, prev), shaped);
            };
    };

// Parameters
maxHold = 0.05;
maxSR = 48000;
att = hslider("att[scale:log]", 0.005*1000, 0.046, maxHold*1000, 0.001)/1000;
att_samples = att*ma.SR:max(1);
half_att_samples = (0.5*att_samples):max(1);
brake_samples = 1.5*att*ma.SR:max(1);
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
