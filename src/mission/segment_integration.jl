"""Sixth-order three-stage Gauss integration of log weight and time.
Implicit stage weights are converged before accepting the interval; the final
endpoint is then evaluated at the integrated weight. Higher order permits long
smooth cruise panels without repeatedly solving dozens of engine states.
"""
function integrate_mission_segment!(rates!, s0, s1, w0, t0;
        rtol=1e-7, max_step=500*nmi_to_m, maxiter=12, validate! = ()->nothing)
    0 < rtol < 1 && max_step > 0 && maxiter > 0 || throw(ArgumentError("Invalid integration controls"))
    s1 >= s0 || throw(ArgumentError("Negative segment distance"))
    w0 > 0 && isfinite(w0) || throw(ArgumentError("Invalid segment weight"))
    f0,v0=rates!(s0,w0)
    isfinite(f0) && f0>=0 && isfinite(v0) && v0>0 || throw(ArgumentError("Invalid segment rates"))
    s1==s0 && return (weight=w0,time=t0,residual=0.,iterations=0)
    panels=max(1,ceil(Int,(s1-s0)/max_step));ds=(s1-s0)/panels
    w=w0;t=t0;worst=0.;iterations=0
    q=sqrt(15.);c1=.5-q/10;c2=.5;c3=.5+q/10
    a11=5/36;a12=2/9-q/15;a13=5/36-q/30
    a21=5/36+q/24;a22=2/9;a23=5/36-q/24
    a31=5/36+q/30;a32=2/9+q/15;a33=5/36
    for panel in 1:panels
        left=s0+(panel-1)*ds;right=s0+panel*ds
        w1=w*exp(-ds*c1*f0);w2=w*exp(-ds*c2*f0);w3=w*exp(-ds*c3*f0)
        converged=false
        for iteration in 1:maxiter
            f1,v1=rates!(left+ds*c1,w1)
            f2,v2=rates!(left+ds*c2,w2)
            f3,v3=rates!(left+ds*c3,w3)
            all(isfinite,(f1,f2,f3,v1,v2,v3)) && min(f1,f2,f3)>=0 && min(v1,v2,v3)>0 ||
                throw(ArgumentError("Invalid segment rates"))
            n1=w*exp(-ds*(a11*f1+a12*f2+a13*f3))
            n2=w*exp(-ds*(a21*f1+a22*f2+a23*f3))
            n3=w*exp(-ds*(a31*f1+a32*f2+a33*f3))
            residual=max(abs(n1/w1-1),abs(n2/w2-1),abs(n3/w3-1));iterations+=1
            if residual<=rtol
                w*=exp(-ds*(5*f1+8*f2+5*f3)/18)
                t+=ds*(5*v1+8*v2+5*v3)/18
                f0,v0=rates!(right,w)
                validate!()
                worst=max(worst,residual);converged=true;break
            end
            w1=n1;w2=n2;w3=n3
        end
        converged || throw(ArgumentError("Segment weight/rate iteration did not converge"))
    end
    (;weight=w,time=t,residual=worst,iterations)
end

# Cache lifetime is one mission iteration: no reuse across aircraft, weight
# iterations or changed integration controls. Entries include every scalar
# trajectory input and the full endpoint state, not just a fuel number.
struct MissionSegmentCacheEntry
    result::NamedTuple{(:weight,:time,:residual,:iterations),Tuple{Float64,Float64,Float64,Int}}
    aero::Vector{Float64}
    engine::Vector{Float64}
    rates::NTuple{3,Float64}
    maximum_temperature::Float64
end
