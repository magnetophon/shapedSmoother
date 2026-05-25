declare name "shapedSmoother";
declare version "0.6";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2025 - 2026, Bart Brouns";
import("stdfaust.lib");

process = test2@att_samples, shapedSmoother(test2);

// ---- LUT configuration ---------------------------------------------------
// 2D table over (slider_value, phase). Indexing by the RAW slider value
// (pre-shapeMap) gives uniform resolution in user-control space. Since
// shapeMap compresses ~99% of the c-range into the top 1% of the slider,
// indexing by c directly causes catastrophic interpolation error at high
// shape — the curves at c=0.984 and c=0.999 are wildly different and a
// linear interp between them at c=0.997 is off by ~50%.
//
// 64 slider bins × 256 phase bins × 4 bytes = 64 KB → L1.
SHAPE_BINS = 64;
PHASE_BINS = 256;

// Function tabulated in slider space; shapeMap evaluated at table-build time.
cheapCurveReleaseFromSlider(s, x) = cheapCurveRelease(shapeMap(s), x);

// 2D lookup, indexed by raw slider value ∈ [0,1] and phase ∈ [0,1].
// C=1 clamps inputs to the table range (cheap insurance at the edges).
curveLUT(s, x) = ba.tabulateNd(1, cheapCurveReleaseFromSlider,
    (SHAPE_BINS, PHASE_BINS,
     0.0, 0.0,
     1.0, 1.0,
     s,   x)).lin;

shapedSmoother(x) = (lookaheadX:env~(_, _, _))
    with {
        lookaheadX = x:ba.slidingMin(att_samples+1, 1+maxHold*maxSR);

        // ---- Control-rate per-shape constants -----------------------------
        // Keep both the raw slider value (for the LUT) and the mapped shape
        // (for everything else — derivativeBase, invPart, etc still need
        // the actual c value).
        attackShapeRaw  = hslider("attack shape",  0, 0, 1, 0.001);
        releaseShapeRaw = hslider("release shape", 0, 0, 1, 0.001);
        attackShape     = shapeMap(attackShapeRaw);
        releaseShape    = shapeMap(releaseShapeRaw);

        attackInvCurveScale  = 1/curveScale(attackShape);
        releaseInvCurveScale = 1/curveScale(releaseShape);

        attackMaxDeriv  = maxDerivativeBaseAttack(attackShape);
        releaseMaxDeriv = maxDerivativeBaseAttack(releaseShape);

        attNRSteps = (att*ma.SR):max(1);
        relNRSteps = (rel*ma.SR):max(1);
        attStep    = 1/attNRSteps;
        relStep    = 1/relNRSteps;

        env(prev, prevPhase, prevTotalStep, lookaheadX) = result, newPhase, totalStep
            with {
                attacking = (prev>lookaheadX)&(att>0);
                releasing = (prev<lookaheadX)&(rel>0);
                active    = attacking|releasing;

                shape         = select2(releasing, attackShape,         releaseShape);
                shapeRaw      = select2(releasing, attackShapeRaw,      releaseShapeRaw);
                invCurveScale = select2(releasing, attackInvCurveScale, releaseInvCurveScale);
                maxDerivVal   = select2(releasing, attackMaxDeriv,      releaseMaxDeriv);
                step          = select2(releasing, attStep,             relStep);

                todo = lookaheadX-prev;
                totalStep = select2(releasing,
                    todo:min(prevTotalStep),
                    todo:max(prevTotalStep)) * active;

                // derivativeBaseAttack(c, x) == derivativeBase(c, 1-x)
                prevY     = select2(releasing, 1-prevPhase, prevPhase);
                prevSpeed = prevTotalStep * derivativeBase(shape, prevY);

                speedRatio   = prevSpeed/totalStep;
                clampedRatio = max(0, min(speedRatio, maxDerivVal));

                cD1     = shape*clampedRatio + 1;
                invPart = sqrt(max(0, 1 - 4*clampedRatio*(1-shape)*cD1));
                invDen2 = 1/(2*cD1);
                topRelease    = (1 + invPart)*invDen2;
                bottomRelease = (1 - invPart)*invDen2;

                // Projected target via LUT, indexed by RAW slider value.
                //   1 - cheapCurveAttack(c, 1-b) == cheapCurveRelease(c, b)
                gonnaDoPhase  = select2(releasing, bottomRelease, topRelease);
                gonnaDoCurve  = curveLUT(shapeRaw, gonnaDoPhase);
                gonnaDoFactor = select2(releasing, gonnaDoCurve, 1 - gonnaDoCurve);
                projected     = gonnaDoFactor*totalStep + prev;
                gonnaMakeIt   = (projected>lookaheadX);

                inner = select2(gonnaMakeIt, bottomRelease, topRelease);
                phaseAtMatchingSpeed = select2(releasing, 1 - inner, inner);

                newPhase = (phaseAtMatchingSpeed+step):min(1-step):max(step)*active;

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
// nuchi's curves, with attack/release as horizontal flips:
//   https://www.wolframalpha.com/input?i=antiderivative+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29
//   https://www.wolframalpha.com/input?i=inverse+function+of+x%281-x%29%2F%28c+*+x%5E2+%2B+1+-+c%29+%2F+C
//   https://www.desmos.com/calculator/9wtfhymvr0
//   scaled c:        https://www.desmos.com/calculator/dynyjjkuli
//   horizontal flip: https://www.desmos.com/calculator/bsr8cdn21v

cheapCurveBase(c, x) = (log(c*(x*x-1)+1)+2*sqrt((1/c)-1)*atan(x/sqrt((1/c)-1))-2*x)/(2*c);
curveScale(c) = cheapCurveBase(c, 1)-cheapCurveBase(c, 0);

cheapCurveRelease(c, x) = (cheapCurveBase(c, x)-cheapCurveBase(c, 0))/curveScale(c);
cheapCurveAttack(c, x)  = 1 - cheapCurveRelease(c, 1-x);

derivativeBase(c, y) = y*(1-y)/(c*y*y + 1 - c);

derivativeBaseRelease(c, x) = derivativeBase(c, x);
derivativeBaseAttack(c, x)  = derivativeBase(c, 1-x);

derivativeRelease(c, x) = derivativeBaseRelease(c, x)/curveScale(c);
derivativeAttack(c, x)  = derivativeBaseAttack(c, x)/curveScale(c);

peakPhaseAttack(c) = 1-sqrt(1-c)*(1-sqrt(1-c))/max(1e-10, c);
maxDerivativeBaseAttack(c) = derivativeBaseAttack(c, peakPhaseAttack(c));
maxDerivativeAttack(c)     = derivativeAttack(c, peakPhaseAttack(c));

inverseDerivativePart(c, D)          = sqrt(max(0, 1-4*D*(1-c)*(c*D+1)));
inverseDerivativeTopRelease(c, D)    = (1+inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeBottomRelease(c, D) = (1-inverseDerivativePart(c, D))/(2*(c*D+1));
inverseDerivativeTopAttack(c, D)     = 1-inverseDerivativeTopRelease(c, D);
inverseDerivativeBottomAttack(c, D)  = 1-inverseDerivativeBottomRelease(c, D);

shapeMap(c) = 1-0.9999*exp(-8.2*pow(c, 1.3));
