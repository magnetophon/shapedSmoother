declare name "shapedSmoother";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2025, Bart Brouns";

import("stdfaust.lib");
import("max_c.lib");

/*


* TODO:
* make self-correcting when using fullDif
* integrate attack and release
* give shapes, derivative and invDer as parameters?
* put everything in lookup tables afap


/*/

// working for hybrid:
// (1024,32)

// max_iter = 1024;
// 128 gives the same result as 1024
max_iter = 128;
// 64 doesn't
// max_iter = 64;
// max_iter = 32;
// max_iter = 16;
// max_iter = 8;
// for diagram:
// max_iter = 4;

// hybrid_iter = 8;
hybrid_iter = 16;
// hybrid_iter = 32;
// hybrid_iter = 64;

// inverseTableSize = 1024;
// inverseTableSize = 2048;
inverseTableSize = pow(2,14);
y_time =
  ((ba.time:min(inverseTableSize)/inverseTableSize));




process =
  // ternary_search_max(512,dwsdx,0,1,shape2c(ba.time:min(nrShapes+1)))
  // -max_c_rdtable(ba.time:min(nrShapes+1))

  (
    (
      dwsdx_hybrid_lut_bisection_inverse(nrShapes, y_time, 0)
      // dwsdx_inverse_lut(nrShapes, y_time,0)
      : dwsdx(shape2c(nrShapes))
    )
    - y_time
  );

// <(ma.EPSILON*hslider("prec", 4, 1, 1024, 1))
// :hbargraph("same", 0, 1)

// ,(             invert_dwsdx(c, hslider("y", 0, 0, 2.95, 0.001),checkbox("top")):hbargraph("[3]x N-R", 0, 1)
// : dwsdx(c):hbargraph("[4]new y N-R", 0, 2.95))

// (max_dwsdx(c):vbargraph("max val grid[unit:dB]", 0, 3))
// ,

// (ternary_search_max(ternary_iter,dwsdx,0,1,c):vbargraph("max val ternary[unit:dB]", 0, 1))
// , (max_c_lut(c):vbargraph("max val table[unit:dB]", 0, 1))
// ,
// ((ternary_search_max(512,dwsdx,0,1,c)==ternary_search_max(ternary_iter,dwsdx,0,1,c)):hbargraph("same", 0, 1))

// TODO: find out why this doesn't always match up:
// (abs(max_c_rdtable(shapeSlider)-ternary_search_max(512,dwsdx,0,1,shape2c(shapeSlider)))<(ma.EPSILON*hslider("prec", 4, 1, 4096, 1)))
// (abs(max_c_rdtable(shapeSlider)-ternary_search_max(512,dwsdx,0,1,shape2c(shapeSlider)))<(0.001*hslider("prec", ma.EPSILON, 0, 1, 0.001)))
// max_c_rdtable(shapeSlider)-ternary_search_max(512,dwsdx,0,1,shape2c(shapeSlider))
// :hbargraph("same", 0, 1)

// (max_c_rdtable(shapeSlider):hbargraph("lut", 0.5, 1))
// ,(ternary_search_max(512,dwsdx,0,1,c):hbargraph("search", 0.5, 1))



// (abs(ternary_search_max(512,dwsdx,0,1,c)-ternary_search_max(ternary_iter,dwsdx,0,1,c))<(ma.EPSILON*hslider("prec", 4, 1, 1024, 1))):hbargraph("same", 0, 1))


// , ((min(ternary_search_max(ternary_iter,dwsdx,0,1,c),max_c_lut(c))
// / max(ternary_search_max(ternary_iter,dwsdx,0,1,c),max_c_lut(c)))
// :hbargraph("error", 0.999, 1.001)
// )
// ;

// testSig,
// shapedSmoother(testSig);

// test2@att_samples,
// shapedSmoother(test2:ba.slidingMax(att_samples,maxHold*maxSR));

// shapedSmoother(x);

// y_target = hslider("y", 0, 0, 1, 0.001)*max_c_lut(c);
// y_target = hslider("y", 0, 0, 1, 0.001)*max_c_rdtable(shapeSlider):si.smoo;
// y_target = hslider("y", 0, 0, 1, 0.001)*max_c_rdtable(shapeSlider);


// c =
// hslider("c", 0.001, 0.001, 3, 0.001)
// * -1
// ;
// shape2c(shapeSlider);
shape2c(slider) =
  slider
  /nrShapes
  <:pow(1+0.42*_)
    * maxShape
    // + ma.EPSILON
    + 0.0000000000001
    // + 0.00000000000001
;
// shapeSlider = hslider("shape", 0, 0, nrShapes, 1):int;
// shapeSlider = hslider("shape", nrShapes, 0, nrShapes, 1):floor;
// TODO: make this a pow of 2 and compensate elsewhere
nrShapes = 127;
// nrShapes = 3;
maxShape = 2.42;
// maxShape = 3;

inverseTargetDerivativeBottom(y) =
  // NR_dwsdx_hybrid_lut_bisection_inverse(0, y, 0);
  dwsdx_hybrid_lut_bisection_inverse(0, y, 0);
// dwsdx_inverse_lut(0, y,0);
// ((0-acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI) ;

inverseTargetDerivativeTop(y) =
  // NR_dwsdx_hybrid_lut_bisection_inverse(0, y, 1);
  dwsdx_hybrid_lut_bisection_inverse(0, y, 1);
// dwsdx_inverse_lut(0, y,1);
// (acos(y/(0.5*ma.PI))+0.5*ma.PI)/ma.PI ;

targetDerivative(x) =
  dwsdx(shape2c(0), x);
// 0.5*ma.PI * cos(ma.PI*x-(0.5*ma.PI));

targetCurve(x) =
  ws(shape2c(0),x);
// (sin((x-0.5)*ma.PI)+1)*0.5;

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

  // which mult factor do we need to match the speed we had previously?
  // plug that into inverseTargetDerivative to get the phase



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


// warp
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


dwsdx_bisection_inverse(c_int, y_target,bottom_top) =
  select2(bottom_top
          // , bisection_inverse_bottom(dwsdx(c), y_target, 0, ternary_search_max(ternary_iter,dwsdx,0,1,c))
          // , bisection_inverse_top(dwsdx(c), y_target, ternary_search_max(ternary_iter,dwsdx,0,1,c),1)
          // , bisection_inverse_bottom(dwsdx(c), y_target, 0,            max_c_lut(c))
          // , bisection_inverse_top   (dwsdx(c), y_target, max_c_lut(c), 1)

          // , bisection_inverse_bottom(dwsdx(c), y_target, 0, max_c_rdtable(shapeSlider))
          // , bisection_inverse_top   (dwsdx(c), y_target, max_c_rdtable(shapeSlider), 1)

         , dwsdx_bisection_inverse_bottom(c_int,y_target)
         , dwsdx_bisection_inverse_top(c_int,y_target)
         )
;
dwsdx_bisection_inverse_bottom(c_int,y_target) =
  bisection_inverse_bottom(max_iter,dwsdx(shape2c(c_int)), y_target, 0, max_c_rdtable(c_int));
dwsdx_bisection_inverse_top(c_int,y_target) =
  bisection_inverse_top(max_iter,dwsdx(shape2c(c_int)), y_target, max_c_rdtable(c_int), 1);

// put dwsdx_bisection_inverse(c, y_target) in a 3d table
dwsdx_inverse_lut(c_int, y_target,bottom_top) =
  ba.tabulateNd(1, dwsdx_bisection_inverse,
                (nrShapes, inverseTableSize, 2,
                 0, 0, 0,
                 nrShapes, 1, 1,
                 c_int,y_target, bottom_top)
               ).cub;
old_dwsdx_inverse_lut(c_int, y_target,bottom_top) =
  select2(bottom_top
         , ba.tabulateNd(0, dwsdx_bisection_inverse_bottom,
                         (nrShapes, inverseTableSize,
                          0, 0,
                          nrShapes, 1,
                          c_int,y_target)
                         // TODO: use 2 times .val to bracket
                        ).val

         , ba.tabulateNd(0, dwsdx_bisection_inverse_top,
                         (nrShapes, inverseTableSize,
                          0, 0,
                          nrShapes, 1,
                          c_int,y_target)
                         // TODO: use 2 times .val to bracket
                        ).val
         );

// get val from table above and from table below the y_target we are looking for
// if dwsdx(c,middle) > y_target
NR_dwsdx_hybrid_lut_bisection_inverse(c_int, y_target,bottom_top) =
  invert_dwsdx_NR(dwsdx_inverse_lut(c_int, y_target,bottom_top),c_int, y_target,bottom_top) ;

dwsdx_hybrid_lut_bisection_inverse(c_int, y_target,bottom_top) =
  select2(bottom_top
          // , bisection_inverse_bottom(dwsdx(shape2c(c_int)), y_target, 0, max_c_rdtable(c_int))
          // , bisection_inverse_top(dwsdx(shape2c(c_int)), y_target, max_c_rdtable(c_int), 1)
         , bisection_inverse_bottom(hybrid_iter,dwsdx(shape2c(c_int)), y_target
                                                                       // error: e-07, with big lut: e-08
                                                                       // , dwsdx_inverse_lut(c_int, y_target-(0.1/inverseTableSize),bottom_top)
                                                                       // , dwsdx_inverse_lut(c_int, y_target+(0.1/inverseTableSize),bottom_top)
                                                                       // error: e-08, but with inverseTableSize = pow(2,14), mostly -e10, except the first: e-09
                                                                       // , dwsdx_inverse_lut(c_int, y_target-(0.01/inverseTableSize),bottom_top)
                                                                       // , dwsdx_inverse_lut(c_int, y_target+(0.01/inverseTableSize),bottom_top)
                                                                       // error: e-07, with big lut: e-11

                                                                       // , dwsdx_inverse_lut(c_int, (y_target-(0.0001/inverseTableSize)):max(0) , bottom_top)
                                                                       // , dwsdx_inverse_lut(c_int, (y_target+(0.0001/inverseTableSize)):min(1) , bottom_top)


                                                                       //TODO: optimize the delta
                                                                       // , dwsdx_inverse_lut(c_int, (y_target-(0.0000001/inverseTableSize)), bottom_top)
                                                                       // , dwsdx_inverse_lut(c_int, (y_target+(0.0000001/inverseTableSize)), bottom_top)
                                    , dwsdx_inverse_lut(c_int, (y_target-(0.00000001/inverseTableSize)), bottom_top)
                                    , dwsdx_inverse_lut(c_int, (y_target+(0.00000001/inverseTableSize)), bottom_top)

                                      // error: e-07, with big lut: e-10
                                      // , dwsdx_inverse_lut(c_int, y_target-(0.00001/inverseTableSize),bottom_top)
                                      // , dwsdx_inverse_lut(c_int, y_target+(0.00001/inverseTableSize),bottom_top)
                                   )
         , bisection_inverse_top(hybrid_iter,dwsdx(shape2c(c_int)), y_target
                                 , dwsdx_inverse_lut(c_int, y_target+(0.01/inverseTableSize),bottom_top)
                                 , dwsdx_inverse_lut(c_int, y_target-(0.01/inverseTableSize),bottom_top)
                                )
         )
  // TODO: remove this fix for 1st val
  <:(select2(
        abs(dwsdx(shape2c(c_int),_)-y_target)
        >
        abs(dwsdx(shape2c(c_int),dwsdx_inverse_lut(c_int, y_target , bottom_top))-y_target)
      , _,
      dwsdx_inverse_lut(c_int, y_target , bottom_top))
    )
;

// Bisection search for inverse: finds x such that f(c,x) == y

bisection_inverse_bottom(iter,f, y, start, end) =
  new_left_right(f, y, start, end)
  : seq(i, iter, new_left_right(f, y))
    :>_/2 // At the end, mid is the best estimate
with {
  // Compute new [left, right] bounds
  new_left_right(f, y, left, right) =
    new_left(f, y, left, right), new_right(f, y, left, right)
  with {
  mid = (left + right) / 2;
  val = f(mid);
  // If val < y, search right half, else search left half
  new_left(f, y, left, right) = select2(val > y, mid, left);
  new_right(f, y, left, right) = select2(val > y, right, mid);
};
};
bisection_inverse_top(iter,f, y, start, end) =
  new_left_right(f, y, start, end)
  : seq(i, iter, new_left_right(f, y))
    :>_/2 // At the end, mid is the best estimate
with {
  // Compute new [left, right] bounds
  new_left_right(f, y, left, right) =
    new_left(f, y, left, right), new_right(f, y, left, right)
  with {
  mid = (left + right) / 2;
  val = f(mid);
  // If val < y, search right half, else search left half
  new_left(f, y, left, right) = select2(val < y, mid, left);
  new_right(f, y, left, right) = select2(val < y, right, mid);
};
};
// Find c such that f(c) == target using the bisection method
bad_bisection_inverse(f, target, start, end) =
  new_left_right(f, target, start, end)
  : seq(i, max_iter, new_left_right(f, target))
    :>_/2  // Use midpoint of last interval as best estimate
with {
  new_left_right(f, target, left, right) =
    new_left(f, target, left, right)
   ,new_right(f, target, left, right)
  with {
  mid = (left + right) / 2;

  // Check sign of f(left)-target vs f(mid)-target to decide next interval
  new_left(f, target, left, right) =
    select2((f(left) - target) * (f(mid) - target) < 0
           , left
           , mid
           );
  new_right(f, target, left, right) =
    select2((f(left) - target) * (f(mid) - target) < 0
           , mid
           , right
           );
};
};
// TODO: maybe we can use a lookup table for the first guess in newton's method
// similarly, can we use a table for thge first guess in the lamb compares?
// TODO: use newton on f and do a hybrid

// find the maximum
max_dwsdx(c) =
  0:seq(i, 800, max_finder(c,i))
with {
  max_finder(c,i,prev_x) =
    max(prev_x, dwsdx(c, i/800));
};

// Ternary/golden-section max search:
//
// only works for Unimodal functions: functions that have a single min or max value
//
// Start with interval [a, b] = [0, 1].
// Evaluate dwsdx at two interior points.
// Based on which is greater, discard part of the interval.
// Repeat until the interval is sufficiently small.
//

// TODO: make c a list, so this can be used with any function.
// add to faustlibraries
ternary_search_max(nr_iter,f,start,end,c) =
  new_left_right(f,c,start,end)
  : seq(i, nr_iter, new_left_right(f,c))
    :>_/2
      // : f(c)
      // : max(f(c),f(c))
with {
  new_left_right(f,c,left, right) =
    new_left(f,c,left, right)
   ,new_right(f,c,left, right)
  with {
  m1 = left + (right - left) / 3;
  m2 = right - (right - left) / 3;
  new_left(f,c,left, right) =
    select2(f(c,m1) < f(c,m2)
           , (left)
           , (m1)
           );
  new_right(f,c,left, right) =
    select2(f(c,m1) < f(c,m2)
           , (m2)
           , (right)
           );
};
};

//tabulate(C, FX, S, r0, r1, x).val
max_c_lut(c) = ba.tabulate(0, ternary_search_max(ternary_iter,dwsdx,0,1), nrShapes+1, shape2c(0), shape2c(nrShapes), c).val;
// max_c_lut(c) = ba.tabulate(0, ternary_search_max(ternary_iter,dwsdx,0,1), nrShapes, shape2c(0), shape2c(nrShapes-1), c).val;

// Golden Section Search for maximizing dwsdx
// doesn't work
golden_section_search_max(f, c, start, end) =
  new_left_right(f, c, start, end)
  : seq(i, max_iter, new_left_right(f, c))
    // :>_/2
  : max(f(c),f(c))
with {
  // Golden ratio constant
  gr = (sqrt(5) + 1) / 2;
  resphi = 2 - gr; // 0.618...

  new_left_right(f, c, left, right) =
    new_left(f, c, left, right),
    new_right(f, c, left, right)
  with {
    d = gr * (right-left);
    m1 = right - resphi * (right - left);
    m2 = left + resphi * (right - left);

    new_left(f, c, left, right) =
      select2(f(c, m1) < f(c, m2), m1, left);

    new_right(f, c, left, right) =
      select2(f(c, m1) < f(c, m2), right, m2);
  };
};


//

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

bracket(C, FX, S, r0, r1, x) = environment {
                                 // Maximum index to access
                                 mid = S-1;

                                 // Create the table
                                 wf = r0 + float(rid(ba.time, 1))*(r1-r0)/float(mid) : FX;

                                 // Prepare the 'float' table read index
                                 id = (x-r0)/(r1-r0)*mid;

                                 // Limit the table read index in [0, mid] if C = 1
                                 rid(x, 0) = x;
                                 rid(x, 1) = max(0, min(x, mid));

                                 // Tabulate an unary 'FX' function on a range [r0, r1]
                                 val = y0 with { y0 = rdtable(S, wf, rid(int(id+0.5), C)); };

                                 // Tabulate an unary 'FX' function over the range [r0, r1] with linear interpolation
                                 lin = it.interpolate_linear(d,y0,y1)
                                 with {
                                   x0 = int(id);
                                   x1 = x0+1;
                                   d  = id-x0;
                                   y0 = rdtable(S, wf, rid(x0, C));
                                   y1 = rdtable(S, wf, rid(x1, C));
                                 };

                                 // Tabulate an unary 'FX' function over the range [r0, r1] with cubic interpolation
                                 cub = it.interpolate_cubic(d,y0,y1,y2,y3)
                                 with {
                                   x0 = x1-1;
                                   x1 = int(id);
                                   x2 = x1+1;
                                   x3 = x2+1;
                                   d  = id-x1;
                                   y0 = rdtable(S, wf, rid(x0, C));
                                   y1 = rdtable(S, wf, rid(x1, C));
                                   y2 = rdtable(S, wf, rid(x2, C));
                                   y3 = rdtable(S, wf, rid(x3, C));
                                 };
                               };

declare tabulate author "Stephane Letz";

// 128 is not enough when maxShape = 99
// singleprecision ternary_iter = 256;
// doubleprecision ternary_iter = 256;
// quadprecision ternary_iter = 256;
// fixedpointprecision ternary_iter = 256;

// for maxShape=3 and nrShapes=15:
// all within 1*ma.EPSILON
ternary_iter = 128;
// ternary_iter = 64;
// ternary_iter = 32;
// # Newton-Raphson solver
invert_dwsdx_NR(x0,c_int, y_target,bottom_top) =
  x0:seq(i, hybrid_iter, newton_block)
with {
  // top starts at 1-ma.EPSILON
  // x0 = select2(bottom_top,ma.EPSILON, 1-ma.EPSILON);
  // x0 = select2(bottom_top,0.001, 0.999);
  // start = hslider("start", ma.EPSILON, ma.EPSILON, 1, 0.001);
  // x0 = select2(bottom_top,start, 1-start);
  newton_block(prev_x) = prev_x + (-F / J)
                         :max(0):min(1)

                                 // <: select2(bottom_top
                                 // , max(0):min(max_c_rdtable(c_int))
                                 // , max(max_c_rdtable(c_int),min(1)))

                                 // :max(0):min(hslider("maximum", 3, 0.5, 3, 0.001))
                                 // :max(0):min(ternary_search_max(dwsdx,0,1,shape2c(c_int)))
  with {
  F = dwsdx(shape2c(c_int), prev_x) - y_target;
  J = second_dwsdx(shape2c(c_int), prev_x);
  // J = select2(checkbox("safe"), unsafe_J, safe_J);
  unsafe_J = second_dwsdx(shape2c(c_int), prev_x);
  // Avoid divide-by-zero
  safe_J = (abs(second_dwsdx(shape2c(c_int), prev_x)):max(ma.EPSILON)) * sign(second_dwsdx(shape2c(c_int), prev_x));
  // ba.signum multiplies by 0 if x is 0
  sign(x) = (x>=0)-(x<0);
};
};
