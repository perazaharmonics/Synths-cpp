
 /*
 *  FilterFactory.hpp
 *  ------------------------------------------------------------------
 *  Unified helpers for creating common IIR / FIR filters and envelopes
 *  in the   **dsp**   namespace.  The file now once again exposes the
 *  complete family requested by the user:
 *
 *    • First‑order One‑Pole LP / HP
 *    • Second‑order biquads: Butterworth, Chebyshev‑I, Bessel, Elliptic
 *    • 3‑pole (≈ 3rd‑order) variants returned as an array<BiQuad,2>
 *    • RBJ parametric EQ sections: BP, Notch, Shelves, PeakEQ
 *    • ADSR envelope and basic envelope followers
 *    • 16‑tap half‑band FIR taps for 2× decimation (unchanged)
 *
 *
 *
  */

#pragma once
#include "OnePole.hpp"
#include "BiQuad.hpp"
#include "EnvelopeADSR.hpp"
#include <array>
#include <vector>
#include <cmath>

namespace sig {

/*==============================  INTERNAL  ==============================*/
namespace detail {

template<typename T>
using BTaps = typename BiQuad<T>::Taps;

template<typename T>
inline void norm(BTaps<T>& t, T a0)
{
    t.b0/=a0; t.b1/=a0; t.b2/=a0;
    t.a1/=a0; t.a2/=a0;
}

/*----- RBJ cookbook primitives we implement here -----------------------*/
template<typename T>
BTaps<T> lpButter(double fs,double fc,double Q)
{
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double a=sw/(2*Q), a0=1+a;
    BTaps<T> t;
    t.b0=(1-cw)*0.5; t.b1=1-cw; t.b2=t.b0;
    t.a1=-2*cw;      t.a2=1-a;
    norm(t,static_cast<T>(a0)); return t;
}

template<typename T>
BTaps<T> hpButter(double fs,double fc,double Q)
{
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double a=sw/(2*Q), a0=1+a;
    BTaps<T> t;
    t.b0=(1+cw)*0.5; t.b1=-(1+cw); t.b2=t.b0;
    t.a1=-2*cw;      t.a2=1-a;
    norm(t,static_cast<T>(a0)); return t;
}

template<typename T>
BTaps<T> bpCook(double fs,double fc,double Q)
{
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double a=sw/(2*Q), a0=1+a;
    BTaps<T> t;
    t.b0= sw*0.5; t.b1=0; t.b2=-sw*0.5;
    t.a1=-2*cw;   t.a2=1-a;
    norm(t,static_cast<T>(a0)); return t;
}

template<typename T>
BTaps<T> notchCook(double fs,double fc,double Q)
{
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double a=sw/(2*Q), a0=1+a;
    BTaps<T> t;
    t.b0=1; t.b1=-2*cw; t.b2=1;
    t.a1=-2*cw; t.a2=1-a;
    norm(t,static_cast<T>(a0)); return t;
}

template<typename T>
BTaps<T> peakCook(double fs,double fc,double gainDB,double Q)
{
    const double A = std::pow(10.0, gainDB/40.0);
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double a=sw/(2*Q);
    BTaps<T> t;
    const double a0 = 1+a/A;
    t.b0 = 1 + a*A;
    t.b1 = -2*cw;
    t.b2 = 1 - a*A;
    t.a1 = -2*cw;
    t.a2 = 1 - a/A;
    norm(t,static_cast<T>(a0)); return t;
}

template<typename T>
BTaps<T> shelfCook(double fs,double fc,double gainDB,double S,bool high)
{
    const double A=std::pow(10.0,gainDB/40.0);
    const double w0=2*M_PI*fc/fs, cw=std::cos(w0), sw=std::sin(w0);
    const double beta=std::sqrt(A)/S, twoBetaSw=2*beta*sw;
    BTaps<T> t; double a0;
    if(!high){
        t.b0=A*((A+1)-(A-1)*cw+twoBetaSw);
        t.b1=2*A*((A-1)-(A+1)*cw);
        t.b2=A*((A+1)-(A-1)*cw-twoBetaSw);
        a0   =(A+1)+(A-1)*cw+twoBetaSw;
        t.a1=-2*((A-1)+(A+1)*cw);
        t.a2=   (A+1)+(A-1)*cw-twoBetaSw;
    }else{
        t.b0=A*((A+1)+(A-1)*cw+twoBetaSw);
        t.b1=-2*A*((A-1)+(A+1)*cw);
        t.b2=A*((A+1)+(A-1)*cw-twoBetaSw);
        a0   =(A+1)-(A-1)*cw+twoBetaSw;
        t.a1= 2*((A-1)-(A+1)*cw);
        t.a2=   (A+1)-(A-1)*cw-twoBetaSw;
    }
    norm(t,static_cast<T>(a0)); return t;
}
} // namespace detail

/*===========================  FILTER FACTORY  ==========================*/
template<typename T=float>
class FilterFactory
{
public:
    /*---- One-pole ----------------------------------------------------*/
    static OnePole<T> OnePoleLP(double fs,double fc){
        OnePole<T> h; h.SetConf(OnePole<T>::Conf::Lowpass); h.Prepare(fs,fc); return h; }
    static OnePole<T> OnePoleHP(double fs,double fc){
        OnePole<T> h; h.SetConf(OnePole<T>::Conf::Highpass); h.Prepare(fs,fc); return h; }

    /*---- Butterworth (RBJ) ------------------------------------------*/
    static BiQuad<T> ButterworthLP(double fs,double fc,double Q){
        BiQuad<T> h; h.SetTaps(detail::lpButter<T>(fs,fc,Q)); return h; }
    static BiQuad<T> ButterworthHP(double fs,double fc,double Q){
        BiQuad<T> h; h.SetTaps(detail::hpButter<T>(fs,fc,Q)); return h; }
    static BiQuad<T> ButterWorthLP(double fs,double fc,double Q){ return ButterworthLP(fs,fc,Q); }
    static BiQuad<T> ButterWorthHP(double fs,double fc,double Q){ return ButterworthHP(fs,fc,Q); }

    /*---- Chebyshev / Bessel / Elliptic (externally implemented) -----*/
    static BiQuad<T> ChebyshevLP (double fs,double fc,double r=0.1,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(sig::ChebyshevLP(fs,fc,r,Q)); return h; }
    static BiQuad<T> ChebyshevHP (double fs,double fc,double r=0.1,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(sig::ChebyshevHP(fs,fc,r,Q)); return h; }

    static BiQuad<T> BesselFirstOrder(double fs,double fc){
        BiQuad<T> h; h.SetTaps(sig::BesselFirstOrder(fs,fc)); return h; }
    static BiQuad<T> Bessel(double fs,double fc,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(sig::Bessel(fs,fc,Q)); return h; }
    static std::array<BiQuad<T>,2> BesselThirdOrder(double fs,double fc,double Q=0.7071){
        auto t=sig::BesselThirdOrder(fs,fc,Q); std::array<BiQuad<T>,2> h;
        h[0].SetTaps(t[0]); h[1].SetTaps(t[1]); return h; }

    static BiQuad<T> EllipticFirstOrder(double fs,double fc){
        BiQuad<T> h; h.SetTaps(sig::EllipticFirstOrder(fs,fc)); return h; }
    static BiQuad<T> Elliptic(double fs,double fc,double r=0.1,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(sig::Elliptic(fs,fc,r,Q)); return h; }
    static std::array<BiQuad<T>,2> EllipticThirdOrder(double fs,double fc,double r=0.1,double Q=0.7071){
        auto t=sig::EllipticThirdOrder(fs,fc,r,Q); std::array<BiQuad<T>,2> h;
        h[0].SetTaps(t[0]); h[1].SetTaps(t[1]); return h; }

    /*---- Parametric EQ sections -------------------------------------*/
    static BiQuad<T> BandPass(double fs,double fc,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(detail::bpCook<T>(fs,fc,Q)); return h; }
    static BiQuad<T> Notch(double fs,double fc,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(detail::notchCook<T>(fs,fc,Q)); return h; }
    static BiQuad<T> PeakEQ(double fs,double fc,double gainDB,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(detail::peakCook<T>(fs,fc,gainDB,Q)); return h; }
    static BiQuad<T> AllPass(double fs,double fc,double Q=0.7071){
        BiQuad<T> h; h.SetTaps(sig::AllPass(fs,fc,Q)); return h; }

    static BiQuad<T> HighShelf(double fs,double fc,double gainDB,double S=1.0){
        BiQuad<T> h; h.SetTaps(detail::shelfCook<T>(fs,fc,gainDB,S,true)); return h; }
    static BiQuad<T> LowShelf (double fs,double fc,double gainDB,double S=1.0){
        BiQuad<T> h; h.SetTaps(detail::shelfCook<T>(fs,fc,gainDB,S,false)); return h; }

    /*---- 10-pole Peak EQ (cascade of 5 biquads) ---------------------*/
    static std::array<BiQuad<T>,5> PeakEQ10Pole(double fs,double fc,double gainDB,double Q=0.7071)
    {
        sig::detail::BTaps<T> t = detail::peakCook<T>(fs,fc,gainDB,Q);
        std::array<BiQuad<T>,5> h;
        for(auto& bq : h) bq.SetTaps(t);
        return h;
    }

    /*---- Envelope utilities ----------------------------------------*/
    static EnvelopeADSR<T> ADSR(double a,double d,double s,double r,double fs){
        EnvelopeADSR<T> env; env.SetAttack(a); env.SetDecay(d);
        env.SetSustain(s); env.SetRelease(r); env.SetSampleRate(fs); return env; }

    static OnePole<T> EnvelopeAttackDetector(double fs,double attMs){
        OnePole<T> h; h.SetConf(OnePole<T>::Conf::Lowpass);
        h.Prepare(fs,1.0/(attMs*0.001*fs)); return h; }

    static OnePole<T> EnvelopeReleaseDetector(double fs,double relMs){
        OnePole<T> h; h.SetConf(OnePole<T>::Conf::Lowpass);
        h.Prepare(fs,1.0/(relMs*0.001*fs)); return h; }

    /*---- FIR 2× half-band decimator taps (poly-phase split) ---------*/
    static std::array<std::vector<T>,2> FIRDecimator(int osFactor=2)
    {
        if(osFactor!=2) return {};
        constexpr T h[16]={-0.0020f,0.0000f, 0.0160f,0.0000f,-0.0800f,0.0000f,
                           0.3000f,0.5000f, 0.3000f,0.0000f,-0.0800f,0.0000f,
                           0.0160f,0.0000f,-0.0020f,0.0000f};
        std::vector<T> even(8),odd(8);
        for(int k=0;k<8;++k){even[k]=h[2*k]; odd[k]=h[2*k+1];}
        return { even, odd };
    }
};

} // namespace sig
