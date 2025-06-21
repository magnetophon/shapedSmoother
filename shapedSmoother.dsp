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
  invert_dwsdx_bottom(2.45, 0.3);

// testSig,
// shapedSmoother(testSig);

// test2@att_samples,
// shapedSmoother(test2:ba.slidingMax(att_samples,maxHold*maxSR));

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

// *************************************** the NEW curves: ******************************

// Based on an algorithm by Dario Sanfilippo:
// https://www.desmos.com/calculator/2hxvf9q194
// Adapted by Bart Brouns:
// https://www.desmos.com/calculator/ubmqgogu2s
// simplified:
// https://www.desmos.com/calculator/cog4ujr7cs
// simplified further:
// https://www.desmos.com/calculator/yiwvcjiony

f(k,x) =
  (1-exp(k*x))
  / (1-exp(k)) ;

fm1m2(c,x) =
  f(-2.42*(c),x);
s(p) = sin(p*ma.PI+1.5*ma.PI)*0.5+0.5;
c2(c,x) = s(fm1m2(c,x));
CurveFormula(c,x) =
  select2(c==0
          // , c2(c,x)
         , c2(c:pow(1+0.42*c),x)
         , s(x)
         )
  :max(0):min(1)
;
Curve(c,x) =
  // CurveFormula(c,x);
  ba.tabulateNd(1, CurveFormula,(nrShapes, tableSize(maxSampleRate),0, 0,1, 1, c,x)).lin;

maxSampleRate = 48000;
// maxSampleRate = 100; // for svg

tableSize(48000) = 1<<16;
tableSize(96000) = 1<<17;
tableSize(sr) = 1<<18;

newCurve(releasing,c,x) =
  select2(releasing
         , Curve(c,x *-1+1 )
           *-1+1
         , Curve(c,x)
         );




// fm1m2(c,x) =
//   f(-2.42*(c),x);
// dfdx(k, x)   = -k * exp(k * x) / (1 - exp(k));
// dfm1m2dx(c, x) = dfdx(-2.42 * c, x);
// f(k, x)      = (1 - exp(k * x)) / (1 - exp(k));
// s(p)         = sin(p * ma.PI + 1.5 * ma.PI) * 0.5 + 0.5;
// dsdp(p)      = 0.5 * ma.PI * cos(p * ma.PI + 1.5 * ma.PI);

// // Derivative of c2 with respect to x
// dc2dx(c, x) = dsdp(fm1m2(c, x)) * dfm1m2dx(c, x);



// simplfied, c a constant number:

// https://www.desmos.com/calculator/icpnkuebe9


f(k, x)      = (1 - exp(k * x)) / (1 - exp(k));
dfdx(k, x)   = -k * exp(k * x) / (1 - exp(k));
second_dfdx(k, x)   = pow(-k,2) * exp(k * x) / (1 - exp(k));
// sine:
s(p)         = sin(p * ma.PI + 1.5 * ma.PI) * 0.5 + 0.5;
dsdp(p)      = 0.5 * ma.PI * cos(p * ma.PI + 1.5 * ma.PI);
second_dsdp(p)= -0.5 * pow(ma.PI,2) * sin(p * ma.PI + 1.5 * ma.PI);
// warped sine:
ws(c,x)        = s(f(c,x));

// Derivative of ws with respect to x
dwsdx(c, x) = dsdp(f(c, x)) * dfdx(c, x);

// second derivative, for use in newton's method:
second_dwsdx(c, x) =
  second_dsdp(f(c, x)) * pow(dfdx(c, x), 2)
  + dsdp(f(c, x)) * second_dfdx(c, x);

max_iter = 3;
// # Newton-Raphson solver
invert_dwsdx_bottom(c, y_target) =
  x0:seq(i, max_iter, newton_block)
with {
  // top starts at 1-ma.EPSILON
  x0 = ma.EPSILON;
  newton_block(prev_x) = prev_x + (-F / J)
  with {
    F = dwsdx(c, prev_x) - y_target;
    J = second_dwsdx(c, prev_x);
  };
};


// # Newton-Raphson solver
// def invert_dwsdx(c, y_target, x0=0.5, max_iter=100, tol=1e-8):
// x = x0
//     for i in range(max_iter):
//     F = dwsdx(c, x) - y_target
//         J = second_dwsdx(c, x)
//             dx = -F / J
//                  x += dx
//                  if abs(dx) < tol:
//                  print(f"Converged after {i+1} iterations.")
//                  break
//                  return x


// TODO: maybe we can use a lookup table for the first guess in newton's method
// similarly, can we use a table for thge first guess in the lamb compares?
// TODO: use newton on f and do a hybrid


/*

// Helper constants and functions
a(c)        = 2.42 * c;
expma(c,x)  = exp(-a(c) * x);
denom(c)    = 1 - exp(-a(c));

// The function itself
c2(c, x) = sin((1 - expma(c,x)) / denom(c) * ma.PI + 1.5 * ma.PI) * 0.5 + 0.5;

// The explicit derivative with respect to x
dc2dx(c, x) =
  0.5 * ma.PI *
  cos(((1 - expma(c,x)) / denom(c)) * ma.PI + 1.5 * ma.PI)
  * a(c) * expma(c,x) / denom(c);

*/
