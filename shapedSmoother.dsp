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
* remove guard-rails in inverseTargetDerivative
* put everything in lookup tables afap


/*/


process(x) =
  // shapedSmoother(test2:ba.slidingMax(att_hold_samples,maxHold*maxSR));
  shapedSmoother(testSig);
// shapedSmoother(x);



shapedSmoother(x) =
  x
 ,(x:att_env~(_,_))
  // , (x:crossfade~(_,_))
  // :(_,!)
  // :(_,targetCurve,targetCurve)
  // f~_
with {
  // crossfade(prev,prevPhase,x) =
  // it.interpolate_linear(naivePhase:targetCurve,prevHold,x)
  // , naivePhase
  att_env(prev,prevPhase,x) =
    select2(attacking,x,
            (speed /totalNRSteps)+prev)
    :min(x)
  , newPhase
    // , gonnaDo(newPhase)
  , fullDif
    // ,(((naivePhase/max(newPhase, ma.EPSILON))-1)*1*attacking)

  with {
  targetCurve(x)= (sin((x-0.5)*ma.PI)+1)*0.5;
  targetDerivative(x) =
    0.5*ma.PI * cos(ma.PI*x-(0.5*ma.PI));

  // which mult factor do we need to match the speed we had previously?
  // plug that into inverseTargetDerivative to get the phase

  // TODO: these guard-rails should not be needed:
  // if we hold the peak they are not needed
  // if we use fullDif they might also not be needed
  inverseTargetDerivativeBottom(y) =
    select2(y>0.5*ma.PI,
            ((0-acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI)
            , 0.5)
  ;

  inverseTargetDerivativeTop(y) =
    select2(y>0.5*ma.PI,
            (acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI
            , 0.5)
  ;


  attacking = (prev<x)*(att>0);
  start_attack = (attacking-attacking'):max(0);
  naivePhase =
    ((_+step)*(1-start_attack):min(1))~_ *attacking;




  // hold = ba.slidingMax(att_hold_samples,maxHold*maxSR,x);
  dif =
    x
    -
    prev
  ;
  fullDif = dif/(1-targetCurve(prevPhase
                               // +(step*hslider("offset", 0.5, -10, 10, 0.001))
                              ));

  // old_totalStep = (x-ba.sAndH(1-attacking,x));
  old_totalStep = (x-prevHold);
  prevHold = ba.sAndH(1-attacking,x);

  totalStep =
    // old_totalStep;
    select2(checkbox("old")
           ,fullDif
           , old_totalStep );


  speed = totalStep* targetDerivative(newPhase):hbargraph("shaper", 0, 2);
  prev_speed =
    totalStep'* targetDerivative(prevPhase);
  //
  //
  //
  // TODO: choose the phase so that at the end of the targetCurve, we are at the val of hold, assuming the phase will keep going up at a contant rate
  // dont match the speed of the prev and add step, but comapare the right things straight away
  // take into account:

  // dTodo = x-prev;
  // dTotal = x-xAtStart;

  // from testing: speed is the same when targetDerivative(prevPhase)*(x-holdMin)' == targetDerivative(newPhase)*(x-holdMin)
  // so targetDerivative(newPhase) == ( targetDerivative(prevPhase)*(x-holdMin)') / (x-holdMin)
  // so newPhase == inverseTargetDerivativeBottom( targetDerivative(prevPhase)*(x-holdMin)') / (x-holdMin)

  // new way, thought-experiment:
  // last_step = prev+ totalStep*step*targetDerivative(newPhase)
  // (x-prev) = prev+ totalStep*step*targetDerivative(newPhase)
  // prev+ totalStep*step*targetDerivative(newPhase) = (x-prev)
  // totalStep*step*targetDerivative(newPhase) = (x-prev) - prev
  // targetDerivative(newPhase) = ((x-prev) - prev)/(totalStep*step)
  // newPhase = inverseTargetDerivativeBottom(((x-prev) - prev)/(totalStep*step))


  // gainCurve = (sin((x - 0.5) * pi) + 1) * 0.5
  // inverseGainCurve(y) = 0.5 + (asin(2*y-1)/pi)
  //
  // gTodo =
  //
  // NEW: match previous speed and then choose a number to add so that gTodo equals gExpected
  // this will differentiate between the two solutions and it should self-correct so we end up at exactly x
  // use the last step as a thought experiment to find that number
  //
  // once we have the formula for the last step:
  // - get the two phases that match the prev speed
  // - use that to pick between the two solutions: calculate the end gain for each, and see which is closer to the actual x
  // - find out how many full steps we would have left if we use that phase
  // - find out how much we have to do in this step: use (NrFullsteps-1)

  gonnaDo(phase) =
    (1-
     targetCurve(phase)
    )*
    totalStep
    +prev
  ;
  /*


      (1-
       targetCurve(phase)
      )*
      totalStep
      +prev

    ==



 */

  newPhase =
    // naivePhase ;
    truePhase;
  truePhase = 
    (phaseAtMatchingSpeed
     // todo: don't blindly add a full step, but directly calculate the phase we need
     // at the moment we match the speed we had before and add a step, but we need to??
     // +thisStep
     +step
     : max(0)
       * attacking)
    :min(1)
     // (_ <: (_*(_<1)))

    :hbargraph("phase", 0, 1)
  ;

  totalNRSteps = att*ma.SR;
  step = 1/(totalNRSteps);
  phaseAtMatchingSpeed =
    select2(
      gonnaMakeIt
    , inverseTargetDerivativeBottom(prev_speed/totalStep)
    , inverseTargetDerivativeTop(prev_speed/totalStep)
    );
  gonnaMakeIt =
    gonnaDo(inverseTargetDerivativeTop(prev_speed/totalStep)
           )>=x
    // )
  ;
  // correctionFactor = (gonnaDo(phaseAtMatchingSpeed
  // +
  // hslider("bool", 0.5, -1, 1, 0.001)*
  // step)/x);
  // get the GRdelta we need to do this sample:
  // - get the whole nr of steps we will do besides this one
  

  // - use that to get the phase
  //   - get the gain we want to be at after this step, so that a whole nr of steps will lead to exactly x
  
  // phaseAfterThis = (prevPhase*totalNRSteps)+1:ceil/totalNRSteps;
  // gonnaDoAfterThis = gonnaDo(phaseAfterThis);
  // deltaXthisStep = prev - gonnaDoAfterThis;
  // deltaXthisStep = x - gonnaDoAfterThis:max(0);
  multFactor = deltaXthisStep * totalStep :hbargraph("mult", 0, 0.5*ma.PI);

  //   - get the mult factor we need to use to end up at that gain
  //   - find the phase using that
  //

  thisStep =
    inverseTargetDerivativeBottom(multFactor);

  // select2(
  // gonnaMakeIt
  // , inverseTargetDerivativeBottom(multFactor)
  // , inverseTargetDerivativeTop(multFactor));
};

  // as long as we're attacking, every time we get anew max, we hold it
  // that's what we use to fade towards
  // if (x>x') && (x>prev) then newPhase
  // else phase +=(1/(att*SR)


};
maxHold = 2;
maxSR = 192000;

att = hslider("att[scale:log]", 0.52, 0.0001, maxHold, 0.001)-0.0001;
att_hold = att*hslider("att hold factor", 1, 0, 1, 0.001) ;
att_hold_samples = att_hold *ma.SR:max(1);

test2 =
  (loop~_)
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(0.75)*5+1);
};
testSig = os.lf_sawpos(0.5)<:(
  ((_>0.25)*hslider("step1", 0.75, -1, 1, 0.001))
  +((_>0.5)*hslider("step2", 0.125, -1, 1, 0.001)));
