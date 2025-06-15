declare name "shapedSmoother";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2025, Bart Brouns";

import("stdfaust.lib");

/*


* TODO:
* make self-correcting when using fullDif
* integrate attack and release
* give shapes, derivative and invDer as parameters?
* put everything in lookup tables afap


/*/


process =
  // testSig,
  // shapedSmoother(testSig);
  test2@att_samples,
  shapedSmoother(test2:ba.slidingMax(att_samples,maxHold*maxSR));
// shapedSmoother(x);



shapedSmoother(x) =
  (x:att_env~(_,_,_))
  :(_,!,!)
with {
  att_env(prev,prevPhase,prevTotalStep,x) =
    select2(attacking,x,
            (speed *step)+prev)
    :min(x)
  , newPhase
  , totalStep
  with {
  targetCurve(x)= (sin((x-0.5)*ma.PI)+1)*0.5;
  targetDerivative(x) =
    0.5*ma.PI * cos(ma.PI*x-(0.5*ma.PI));

  // which mult factor do we need to match the speed we had previously?
  // plug that into inverseTargetDerivative to get the phase

  inverseTargetDerivativeBottom(y) =
    ((0-acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI) ;

  inverseTargetDerivativeTop(y) =
    (acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI ;


  attacking = (prev<x)*(att>0);

  totalStep =
    (x-prev)
    : max(prevTotalStep)
      * attacking;


  speed = totalStep* targetDerivative(newPhase)
          // :hbargraph("shaper", 0, 2)
  ;
  prevSpeed =
    prevTotalStep * targetDerivative(prevPhase+step);

  gonnaDo(phase) =
    (1-
     targetCurve(phase)
    )*
    totalStep
  ;

  newPhase =
    (phaseAtMatchingSpeed
     : max(0)
       * attacking
    )
    :min(1)
     // :hbargraph("phase", 0, 1)
  with {
    phaseAtMatchingSpeed =
      select2(
        gonnaMakeIt
      , inverseTargetDerivativeBottom(prevSpeed/totalStep)
      , inverseTargetDerivativeTop(prevSpeed/totalStep)
      );
    gonnaMakeIt =
      gonnaDo(inverseTargetDerivativeTop(prevSpeed/totalStep))
      +prev

      >=x
      // )
    ;
  };


  totalNRSteps = att*ma.SR;
  step = (1/(totalNRSteps));
};

};
maxHold = 2;
maxSR = 192000;
// maxSR = 48000;

att = hslider("att[scale:log]", 0.2, 0.0001, maxHold, 0.001)-0.0001;
att_hold = att*hslider("att hold factor", 1, 0, 1, 0.001) ;
att_hold_samples = att_hold *ma.SR:max(1);
att_samples = att *ma.SR:max(1);

test2 =
  it.interpolate_linear(
    hslider("noise level", 0, 0, 1, 0.001),
    (loop~_)
    ,
      no.lfnoise(hslider("noise rate", 42, 1, 1000, 1))

  )
with {
  loop(prev) = no.lfnoise0(abs(prev*69)%9:pow(0.75)*5+1);
};
testSig = os.lf_sawpos(0.5)<:(
  ((_>0.25)*hslider("step1", 0.75, -1, 1, 0.001))
  +((_>0.5)*hslider("step2", 0.125, -1, 1, 0.001)));
