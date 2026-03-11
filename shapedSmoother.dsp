declare name "shapedSmoother";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2025, Bart Brouns";

import("stdfaust.lib");

process =
  test2@(att_samples),
  shapedSmoother(test2:ba.slidingMin(att_samples+1,1+maxHold*maxSR));

shapedSmoother(x) =
  (x:att_env~(_,_,_))
  :(_,_,!,_)
with {
  att_env(prev,prevPhase,prevTotalStep,x) =
    select2(attacking,x,
            (speed*step)+prev)
    :max(x)
     // :max(test2@att_samples)
  , newPhase
  , totalStep
  , x
  with {
  shape = shapeMap(hslider("shape", 0, 0, 1, 0.001));
  attacking = (prev>x)*(att>0);

    totalStep =
      (x-prev)
      : min(prevTotalStep)
        * attacking;

    safeTotalStep = select2(abs(totalStep) > 1e-30, 1, totalStep);
    speedRatio = prevSpeed / safeTotalStep * (abs(totalStep) > 1e-30);
    clampedRatio = max(0, min(speedRatio, maxDerivativeAttack(shape)));

    speed = totalStep * derivativeAttack(shape,newPhase);
    prevSpeed = prevTotalStep * derivativeAttack(shape,prevPhase);

    gonnaDo(phase) =
      (1 - cheapCurveAttack(shape,phase)) * totalStep;

    phaseAtMatchingSpeed =
      select2(gonnaMakeIt
      , inverseDerivativeBottomAttack(shape,clampedRatio)
      , inverseDerivativeTopAttack(shape,clampedRatio)
      );

    gonnaMakeIt =
      gonnaDo(inverseDerivativeBottomAttack(shape,clampedRatio))
      + prev
      > x;

    newPhase =
      (phaseAtMatchingSpeed + step)
      :min(1-step)
      :max(step)
       * attacking;

    totalNRSteps = att*ma.SR;
    step = 1/totalNRSteps;
  };
};

// Parameters

maxHold = 0.05;
// maxSR = 192000;
maxSR = 48000;

att = hslider("att[scale:log]", 0.005*1000, 1, maxHold*1000, 0.001)/1000;
att_samples = att * ma.SR : max(1);

// Test signal

test2 =
  it.interpolate_linear(
    hslider("noise level", 0, 0, 1, 0.001),
    (loop~_),
    no.lfnoise(hslider("noise rate", 42, 1, 1000, 1))
  )
with {
  loop(prev) = no.lfnoise0(abs(prev*69)%9:pow(0.75)*5+1);
};
testSig = os.lf_sawpos(0.5)<:(
  ((_>0.25)*hslider("step1", 0.75, -1, 1, 0.001))
  +((_>0.5)*hslider("step2", 0.125, -1, 1, 0.001)));

// *************************************** the NEW curves: ******************************

// New curves by nuchi:
// https://www.wolframalpha.com/input?i=antiderivative+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29
// https://www.wolframalpha.com/input?i=inverse+function+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29+%2F+C
// https://www.desmos.com/calculator/9wtfhymvr0
// with scaling for c:
// https://www.desmos.com/calculator/dynyjjkuli
// with flip horizontal:
// https://www.desmos.com/calculator/bsr8cdn21v



cheapCurveBase(c,x) =
  (log(c*(pow(x,2)-1)+1) + 2*sqrt((1/c)-1) * atan(x/sqrt((1/c)-1)) - 2*x) / (2*c);

curveScale(c) = cheapCurveBase(c,1) - cheapCurveBase(c,0);

cheapCurveRelease(c,x) = (cheapCurveBase(c,x) - cheapCurveBase(c,0)) / curveScale(c);
cheapCurveAttack(c,x) = cheapCurveRelease(c, x*-1+1) * -1 + 1;

derivativeBaseRelease(c,x) = x*(1-x) / (c*pow(x,2) + 1 - c);
derivativeBaseAttack(c,x) = x*(1-x) / (c*pow(x,2) + 1 - (2*c*x));

derivativeRelease(c,x) = derivativeBaseRelease(c,x) / curveScale(c);
derivativeAttack(c,x) = derivativeBaseAttack(c,x) / curveScale(c);

peakPhaseAttack(c) = 1 - sqrt(1-c) * (1 - sqrt(1-c)) / max(1e-10, c);
maxDerivativeAttack(c) = derivativeAttack(c, peakPhaseAttack(c));

inverseDerivativePart(c,x) =
  sqrt(max(0, 1 - 4*(curveScale(c)*x - c*curveScale(c)*x) * (c*curveScale(c)*x + 1)));

inverseDerivativeTopRelease(c,x) = (1 + inverseDerivativePart(c,x)) / (2*(c*curveScale(c)*x + 1));
inverseDerivativeBottomRelease(c,x) = (1 - inverseDerivativePart(c,x)) / (2*(c*curveScale(c)*x + 1));
inverseDerivativeTopAttack(c,x) = inverseDerivativeTopRelease(c,x) * -1 + 1;
inverseDerivativeBottomAttack(c,x) = inverseDerivativeBottomRelease(c,x) * -1 + 1;

shapeMap(c) = 1 - 0.9999 * exp(-8.2 * pow(c, 1.3));

newCurve(releasing,c,x) =
  select2(releasing
  , Curve(c, x*-1+1) * -1 + 1
  , Curve(c, x)
  );
