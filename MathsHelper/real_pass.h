/* fft/real_pass.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#include "gsl_complex.h"
#define ATOMIC double
#define VECTOR(a,stride,i) ((a)[(stride)*(i)])

static void fft_real_pass_2 (const double in[],
                                       const size_t istride,
                                       double out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                             const gsl_complex twiddle[]){
    size_t k, k1;
    
    const size_t factor = 2;
    const size_t m = n / factor;
    const size_t q = n / product;
    const size_t product_1 = product / factor;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1;
        const size_t from1 = from0 + m;
        
        const ATOMIC r0 = VECTOR(in,istride,from0);
        const ATOMIC r1 = VECTOR(in,istride,from1);
        
        const ATOMIC s0 = r0 + r1;
        const ATOMIC s1 = r0 - r1;
        
        const size_t to0 = product * k1;
        const size_t to1 = to0 + product - 1;
        
        VECTOR(out,ostride,to0) = s0;
        VECTOR(out,ostride,to1) = s1;
    }
    
    if (product_1 == 1)
        return;
    
    for (k = 1; k < (product_1 + 1) / 2; k++)
    {
        
        /* forward real transform: w -> conjugate(w) */
        const ATOMIC w_real = GSL_REAL(twiddle[k - 1]);
        const ATOMIC w_imag = -GSL_IMAG(twiddle[k - 1]);
        
        for (k1 = 0; k1 < q; k1++)
        {
            const size_t from0 = k1 * product_1 + 2 * k - 1;
            const size_t from1 = from0 + m;
            
            const ATOMIC f0_real = VECTOR(in,istride,from0);
            const ATOMIC f0_imag = VECTOR(in,istride,from0 + 1);
            
            const ATOMIC f1_real = VECTOR(in,istride,from1);
            const ATOMIC f1_imag = VECTOR(in,istride,from1 + 1);
            
            const ATOMIC z0_real = f0_real;
            const ATOMIC z0_imag = f0_imag;
            
            const ATOMIC z1_real = w_real * f1_real - w_imag * f1_imag;
            const ATOMIC z1_imag = w_real * f1_imag + w_imag * f1_real;
            
            /* compute x = W(2) z */
            
            /* x0 = z0 + z1 */
            const ATOMIC x0_real = z0_real + z1_real;
            const ATOMIC x0_imag = z0_imag + z1_imag;
            
            /* x1 = z0 - z1 */
            const ATOMIC x1_real = z0_real - z1_real;
            const ATOMIC x1_imag = z0_imag - z1_imag;
            
            const size_t to0 = k1 * product + 2 * k - 1;
            const size_t to1 = k1 * product + product - 2 * k - 1;
            
            VECTOR(out,ostride,to0) = x0_real;
            VECTOR(out,ostride,to0 + 1) = x0_imag;
            
            /* stored in conjugate location */
            VECTOR(out,ostride,to1) = x1_real;
            VECTOR(out,ostride,to1 + 1) = -x1_imag;
        }
    }
    
    if (product_1 % 2 == 1)
        return;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1 + product_1 - 1;
        const size_t from1 = from0 + m;
        const size_t to0 = k1 * product + product_1 - 1;
        
        VECTOR(out,ostride,to0) = VECTOR(in,istride,from0);
        VECTOR(out,ostride,to0 + 1) = -VECTOR(in,istride,from1);
    }
    return;

}

static void fft_real_pass_3 (const double in[],
                                       const size_t istride,
                                       double out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                             const gsl_complex twiddle2[]){
    size_t k, k1;
    
    const size_t factor = 3;
    const size_t m = n / factor;
    const size_t q = n / product;
    const size_t product_1 = product / factor;
    
    const ATOMIC tau = sqrt (3.0) / 2.0;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        
        const ATOMIC z0_real = VECTOR(in,istride,from0);
        const ATOMIC z1_real = VECTOR(in,istride,from1);
        const ATOMIC z2_real = VECTOR(in,istride,from2);
        
        const ATOMIC t1 = z1_real + z2_real;
        
        const ATOMIC x0_real = z0_real + t1;
        const ATOMIC x1_real = z0_real - t1 / 2.0;
        const ATOMIC x1_imag = -tau * (z1_real - z2_real);
        
        const size_t to0 = product * k1;
        const size_t to1 = to0 + 2 * product_1 - 1;
        
        VECTOR(out,ostride,to0) = x0_real;
        VECTOR(out,ostride,to1) = x1_real;
        VECTOR(out,ostride,to1 + 1) = x1_imag;
    }
    
    if (product_1 == 1)
        return;
    
    for (k = 1; k < (product_1 + 1) / 2; k++)
    {
        const ATOMIC w1_real = GSL_REAL(twiddle1[k - 1]);
        const ATOMIC w1_imag = -GSL_IMAG(twiddle1[k - 1]);
        const ATOMIC w2_real = GSL_REAL(twiddle2[k - 1]);
        const ATOMIC w2_imag = -GSL_IMAG(twiddle2[k - 1]);
        
        for (k1 = 0; k1 < q; k1++)
        {
            const size_t from0 = k1 * product_1 + 2 * k - 1;
            const size_t from1 = from0 + m;
            const size_t from2 = from1 + m;
            
            const ATOMIC f0_real = VECTOR(in,istride,from0);
            const ATOMIC f0_imag = VECTOR(in,istride,from0 + 1);
            const ATOMIC f1_real = VECTOR(in,istride,from1);
            const ATOMIC f1_imag = VECTOR(in,istride,from1 + 1);
            const ATOMIC f2_real = VECTOR(in,istride,from2);
            const ATOMIC f2_imag = VECTOR(in,istride,from2 + 1);
            
            const ATOMIC z0_real = f0_real;
            const ATOMIC z0_imag = f0_imag;
            const ATOMIC z1_real = w1_real * f1_real - w1_imag * f1_imag;
            const ATOMIC z1_imag = w1_real * f1_imag + w1_imag * f1_real;
            const ATOMIC z2_real = w2_real * f2_real - w2_imag * f2_imag;
            const ATOMIC z2_imag = w2_real * f2_imag + w2_imag * f2_real;
            
            /* compute x = W(3) z */
            
            /* t1 = z1 + z2 */
            const ATOMIC t1_real = z1_real + z2_real;
            const ATOMIC t1_imag = z1_imag + z2_imag;
            
            /* t2 = z0 - t1/2 */
            const ATOMIC t2_real = z0_real - t1_real / 2;
            const ATOMIC t2_imag = z0_imag - t1_imag / 2;
            
            /* t3 = (+/-) sin(pi/3)*(z1 - z2) */
            const ATOMIC t3_real = -tau * (z1_real - z2_real);
            const ATOMIC t3_imag = -tau * (z1_imag - z2_imag);
            
            /* x0 = z0 + t1 */
            const ATOMIC x0_real = z0_real + t1_real;
            const ATOMIC x0_imag = z0_imag + t1_imag;
            
            /* x1 = t2 + i t3 */
            const ATOMIC x1_real = t2_real - t3_imag;
            const ATOMIC x1_imag = t2_imag + t3_real;
            
            /* x2 = t2 - i t3 */
            const ATOMIC x2_real = t2_real + t3_imag;
            const ATOMIC x2_imag = t2_imag - t3_real;
            
            /* apply twiddle factors */
            
            const size_t to0 = k1 * product + 2 * k - 1;
            const size_t to1 = to0 + 2 * product_1;
            const size_t to2 = 2 * product_1 - 2 * k + k1 * product - 1;
            
            /* to0 = 1 * x0 */
            VECTOR(out,ostride,to0) = x0_real;
            VECTOR(out,ostride,to0 + 1) = x0_imag;
            
            /* to1 = 1 * x1 */
            VECTOR(out,ostride,to1) = x1_real;
            VECTOR(out,ostride,to1 + 1) = x1_imag;
            
            /* to2 = 1 * x2 */
            VECTOR(out,ostride,to2) = x2_real;
            VECTOR(out,ostride,to2 + 1) = -x2_imag;
        }
    }
    
    if (product_1 % 2 == 1)
        return;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1 + product_1 - 1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        
        const ATOMIC z0_real = VECTOR(in,istride,from0);
        const ATOMIC z1_real = VECTOR(in,istride,from1);
        const ATOMIC z2_real = VECTOR(in,istride,from2);
        
        const ATOMIC t1 = z1_real - z2_real;
        const ATOMIC x0_real = z0_real + t1 / 2.0;
        const ATOMIC x0_imag = -tau * (z1_real + z2_real);
        const ATOMIC x1_real = z0_real - t1;
        
        const size_t to0 = k1 * product + product_1 - 1;
        const size_t to1 = to0 + 2 * product_1;
        
        VECTOR(out,ostride,to0) = x0_real;
        VECTOR(out,ostride,to0 + 1) = x0_imag;
        VECTOR(out,ostride,to1) = x1_real;
    }
    
    return;

}

static void fft_real_pass_4 (const double in[],
                                       const size_t istride,
                                       double out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                                       const gsl_complex twiddle2[],
                             const gsl_complex twiddle3[]){
    size_t k, k1;
    
    const size_t factor = 4;
    const size_t m = n / factor;
    const size_t q = n / product;
    const size_t product_1 = product / factor;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        const size_t from3 = from2 + m;
        
        const ATOMIC z0_real = VECTOR(in,istride,from0);
        const ATOMIC z1_real = VECTOR(in,istride,from1);
        const ATOMIC z2_real = VECTOR(in,istride,from2);
        const ATOMIC z3_real = VECTOR(in,istride,from3);
        
        /* compute x = W(4) z */
        
        /* t1 = z0 + z2 */
        const ATOMIC t1_real = z0_real + z2_real;
        
        /* t2 = z1 + z3 */
        const ATOMIC t2_real = z1_real + z3_real;
        
        /* t3 = z0 - z2 */
        const ATOMIC t3_real = z0_real - z2_real;
        
        /* t4 = - (z1 - z3) */
        const ATOMIC t4_real = -(z1_real - z3_real);
        
        /* x0 = t1 + t2 */
        const ATOMIC x0_real = t1_real + t2_real;
        
        /* x1 = t3 + i t4 */
        const ATOMIC x1_real = t3_real;
        const ATOMIC x1_imag = t4_real;
        
        /* x2 = t1 - t2 */
        const ATOMIC x2_real = t1_real - t2_real;
        
        const size_t to0 = product * k1;
        const size_t to1 = to0 + 2 * product_1 - 1;
        const size_t to2 = to1 + 2 * product_1;
        
        VECTOR(out,ostride,to0) = x0_real;
        VECTOR(out,ostride,to1) = x1_real;
        VECTOR(out,ostride,to1 + 1) = x1_imag;
        VECTOR(out,ostride,to2) = x2_real;
    }
    
    if (product_1 == 1)
        return;
    
    for (k = 1; k < (product_1 + 1) / 2; k++)
    {
        ATOMIC w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;
        w1_real = GSL_REAL(twiddle1[k - 1]);
        w1_imag = -GSL_IMAG(twiddle1[k - 1]);
        w2_real = GSL_REAL(twiddle2[k - 1]);
        w2_imag = -GSL_IMAG(twiddle2[k - 1]);
        w3_real = GSL_REAL(twiddle3[k - 1]);
        w3_imag = -GSL_IMAG(twiddle3[k - 1]);
        
        for (k1 = 0; k1 < q; k1++)
        {
            const size_t from0 = k1 * product_1 + 2 * k - 1;
            const size_t from1 = from0 + m;
            const size_t from2 = from1 + m;
            const size_t from3 = from2 + m;
            
            const ATOMIC f0_real = VECTOR(in,istride,from0);
            const ATOMIC f0_imag = VECTOR(in,istride,from0 + 1);
            const ATOMIC f1_real = VECTOR(in,istride,from1);
            const ATOMIC f1_imag = VECTOR(in,istride,from1 + 1);
            const ATOMIC f2_real = VECTOR(in,istride,from2);
            const ATOMIC f2_imag = VECTOR(in,istride,from2 + 1);
            const ATOMIC f3_real = VECTOR(in,istride,from3);
            const ATOMIC f3_imag = VECTOR(in,istride,from3 + 1);
            
            const ATOMIC z0_real = f0_real;
            const ATOMIC z0_imag = f0_imag;
            const ATOMIC z1_real = w1_real * f1_real - w1_imag * f1_imag;
            const ATOMIC z1_imag = w1_real * f1_imag + w1_imag * f1_real;
            const ATOMIC z2_real = w2_real * f2_real - w2_imag * f2_imag;
            const ATOMIC z2_imag = w2_real * f2_imag + w2_imag * f2_real;
            const ATOMIC z3_real = w3_real * f3_real - w3_imag * f3_imag;
            const ATOMIC z3_imag = w3_real * f3_imag + w3_imag * f3_real;
            
            /* compute x = W(4) z */
            
            /* t1 = z0 + z2 */
            const ATOMIC t1_real = z0_real + z2_real;
            const ATOMIC t1_imag = z0_imag + z2_imag;
            
            /* t2 = z1 + z3 */
            const ATOMIC t2_real = z1_real + z3_real;
            const ATOMIC t2_imag = z1_imag + z3_imag;
            
            /* t3 = z0 - z2 */
            const ATOMIC t3_real = z0_real - z2_real;
            const ATOMIC t3_imag = z0_imag - z2_imag;
            
            /* t4 = - (z1 - z3) */
            const ATOMIC t4_real = -(z1_real - z3_real);
            const ATOMIC t4_imag = -(z1_imag - z3_imag);
            
            /* x0 = t1 + t2 */
            const ATOMIC x0_real = t1_real + t2_real;
            const ATOMIC x0_imag = t1_imag + t2_imag;
            
            /* x1 = t3 + i t4 */
            const ATOMIC x1_real = t3_real - t4_imag;
            const ATOMIC x1_imag = t3_imag + t4_real;
            
            /* x2 = t1 - t2 */
            const ATOMIC x2_real = t1_real - t2_real;
            const ATOMIC x2_imag = t1_imag - t2_imag;
            
            /* x3 = t3 - i t4 */
            const ATOMIC x3_real = t3_real + t4_imag;
            const ATOMIC x3_imag = t3_imag - t4_real;
            
            const size_t to0 = k1 * product + 2 * k - 1;
            const size_t to1 = to0 + 2 * product_1;
            const size_t to2 = 2 * product_1 - 2 * k + k1 * product - 1;
            const size_t to3 = to2 + 2 * product_1;
            
            VECTOR(out,ostride,to0) = x0_real;
            VECTOR(out,ostride,to0 + 1) = x0_imag;
            
            VECTOR(out,ostride,to1) = x1_real;
            VECTOR(out,ostride,to1 + 1) = x1_imag;
            
            VECTOR(out,ostride,to3) = x2_real;
            VECTOR(out,ostride,to3 + 1) = -x2_imag;
            
            VECTOR(out,ostride,to2) = x3_real;
            VECTOR(out,ostride,to2 + 1) = -x3_imag;
        }
    }
    
    if (product_1 % 2 == 1)
        return;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1 + product_1 - 1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        const size_t from3 = from2 + m;
        
        const ATOMIC x0 = VECTOR(in,istride,from0);
        const ATOMIC x1 = VECTOR(in,istride,from1);
        const ATOMIC x2 = VECTOR(in,istride,from2);
        const ATOMIC x3 = VECTOR(in,istride,from3);
        
        const ATOMIC t1 = (1.0 / sqrt (2.0)) * (x1 - x3);
        const ATOMIC t2 = (1.0 / sqrt (2.0)) * (x1 + x3);
        
        const size_t to0 = k1 * product + 2 * k - 1;
        const size_t to1 = to0 + 2 * product_1;
        
        VECTOR(out,ostride,to0) = x0 + t1;
        VECTOR(out,ostride,to0 + 1) = -x2 - t2;
        
        VECTOR(out,ostride,to1) = x0 - t1;
        VECTOR(out,ostride,to1 + 1) = x2 - t2;
    }
    return;
}

static void fft_real_pass_5 (const double in[],
                                       const size_t istride,
                                       double out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                                       const gsl_complex twiddle2[],
                                       const gsl_complex twiddle3[],
                             const gsl_complex twiddle4[]){
    size_t k, k1;
    
    const size_t factor = 5;
    const size_t m = n / factor;
    const size_t q = n / product;
    const size_t product_1 = product / factor;
    
    const ATOMIC sina = sin (2.0 * M_PI / 5.0);
    const ATOMIC sinb = sin (2.0 * M_PI / 10.0);
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        const size_t from3 = from2 + m;
        const size_t from4 = from3 + m;
        
        const ATOMIC z0_real = VECTOR(in,istride,from0);
        const ATOMIC z1_real = VECTOR(in,istride,from1);
        const ATOMIC z2_real = VECTOR(in,istride,from2);
        const ATOMIC z3_real = VECTOR(in,istride,from3);
        const ATOMIC z4_real = VECTOR(in,istride,from4);
        
        /* t1 = z1 + z4 */
        const ATOMIC t1_real = z1_real + z4_real;
        
        /* t2 = z2 + z3 */
        const ATOMIC t2_real = z2_real + z3_real;
        
        /* t3 = z1 - z4 */
        const ATOMIC t3_real = z1_real - z4_real;
        
        /* t4 = z2 - z3 */
        const ATOMIC t4_real = z2_real - z3_real;
        
        /* t5 = t1 + t2 */
        const ATOMIC t5_real = t1_real + t2_real;
        
        /* t6 = (sqrt(5)/4)(t1 - t2) */
        const ATOMIC t6_real = (sqrt (5.0) / 4.0) * (t1_real - t2_real);
        
        /* t7 = z0 - ((t5)/4) */
        const ATOMIC t7_real = z0_real - t5_real / 4.0;
        
        /* t8 = t7 + t6 */
        const ATOMIC t8_real = t7_real + t6_real;
        
        /* t9 = t7 - t6 */
        const ATOMIC t9_real = t7_real - t6_real;
        
        /* t10 = -(sin(2 pi/5) t3 + sin(2 pi/10) t4 ) */
        const ATOMIC t10_real = -sina * t3_real - sinb * t4_real;
        
        /* t11 = -(sin(2 pi/10) t3 - sin(2 pi/5) t4) */
        const ATOMIC t11_real = -sinb * t3_real + sina * t4_real;
        
        /* x0 = z0 + t5 */
        const ATOMIC x0_real = z0_real + t5_real;
        
        /* x1 = t8 + i t10 */
        const ATOMIC x1_real = t8_real;
        const ATOMIC x1_imag = t10_real;
        
        /* x2 = t9 + i t11 */
        const ATOMIC x2_real = t9_real;
        const ATOMIC x2_imag = t11_real;
        
        const size_t to0 = product * k1;
        const size_t to1 = to0 + 2 * product_1 - 1;
        const size_t to2 = to1 + 2 * product_1;
        
        VECTOR(out,ostride,to0) = x0_real;
        VECTOR(out,ostride,to1) = x1_real;
        VECTOR(out,ostride,to1 + 1) = x1_imag;
        VECTOR(out,ostride,to2) = x2_real;
        VECTOR(out,ostride,to2 + 1) = x2_imag;
    }
    
    if (product_1 == 1)
        return;
    
    for (k = 1; k < (product_1 + 1) / 2; k++)
    {
        const ATOMIC w1_real = GSL_REAL(twiddle1[k - 1]);
        const ATOMIC w1_imag = -GSL_IMAG(twiddle1[k - 1]);
        const ATOMIC w2_real = GSL_REAL(twiddle2[k - 1]);
        const ATOMIC w2_imag = -GSL_IMAG(twiddle2[k - 1]);
        const ATOMIC w3_real = GSL_REAL(twiddle3[k - 1]);
        const ATOMIC w3_imag = -GSL_IMAG(twiddle3[k - 1]);
        const ATOMIC w4_real = GSL_REAL(twiddle4[k - 1]);
        const ATOMIC w4_imag = -GSL_IMAG(twiddle4[k - 1]);
        
        for (k1 = 0; k1 < q; k1++)
        {
            const size_t from0 = k1 * product_1 + 2 * k - 1;
            const size_t from1 = from0 + m;
            const size_t from2 = from1 + m;
            const size_t from3 = from2 + m;
            const size_t from4 = from3 + m;
            
            const ATOMIC f0_real = VECTOR(in,istride,from0);
            const ATOMIC f0_imag = VECTOR(in,istride,from0 + 1);
            const ATOMIC f1_real = VECTOR(in,istride,from1);
            const ATOMIC f1_imag = VECTOR(in,istride,from1 + 1);
            const ATOMIC f2_real = VECTOR(in,istride,from2);
            const ATOMIC f2_imag = VECTOR(in,istride,from2 + 1);
            const ATOMIC f3_real = VECTOR(in,istride,from3);
            const ATOMIC f3_imag = VECTOR(in,istride,from3 + 1);
            const ATOMIC f4_real = VECTOR(in,istride,from4);
            const ATOMIC f4_imag = VECTOR(in,istride,from4 + 1);
            
            const ATOMIC z0_real = f0_real;
            const ATOMIC z0_imag = f0_imag;
            const ATOMIC z1_real = w1_real * f1_real - w1_imag * f1_imag;
            const ATOMIC z1_imag = w1_real * f1_imag + w1_imag * f1_real;
            const ATOMIC z2_real = w2_real * f2_real - w2_imag * f2_imag;
            const ATOMIC z2_imag = w2_real * f2_imag + w2_imag * f2_real;
            const ATOMIC z3_real = w3_real * f3_real - w3_imag * f3_imag;
            const ATOMIC z3_imag = w3_real * f3_imag + w3_imag * f3_real;
            const ATOMIC z4_real = w4_real * f4_real - w4_imag * f4_imag;
            const ATOMIC z4_imag = w4_real * f4_imag + w4_imag * f4_real;
            
            /* compute x = W(5) z */
            
            /* t1 = z1 + z4 */
            const ATOMIC t1_real = z1_real + z4_real;
            const ATOMIC t1_imag = z1_imag + z4_imag;
            
            /* t2 = z2 + z3 */
            const ATOMIC t2_real = z2_real + z3_real;
            const ATOMIC t2_imag = z2_imag + z3_imag;
            
            /* t3 = z1 - z4 */
            const ATOMIC t3_real = z1_real - z4_real;
            const ATOMIC t3_imag = z1_imag - z4_imag;
            
            /* t4 = z2 - z3 */
            const ATOMIC t4_real = z2_real - z3_real;
            const ATOMIC t4_imag = z2_imag - z3_imag;
            
            /* t5 = t1 + t2 */
            const ATOMIC t5_real = t1_real + t2_real;
            const ATOMIC t5_imag = t1_imag + t2_imag;
            
            /* t6 = (sqrt(5)/4)(t1 - t2) */
            const ATOMIC t6_real = (sqrt (5.0) / 4.0) * (t1_real - t2_real);
            const ATOMIC t6_imag = (sqrt (5.0) / 4.0) * (t1_imag - t2_imag);
            
            /* t7 = z0 - ((t5)/4) */
            const ATOMIC t7_real = z0_real - t5_real / 4.0;
            const ATOMIC t7_imag = z0_imag - t5_imag / 4.0;
            
            /* t8 = t7 + t6 */
            const ATOMIC t8_real = t7_real + t6_real;
            const ATOMIC t8_imag = t7_imag + t6_imag;
            
            /* t9 = t7 - t6 */
            const ATOMIC t9_real = t7_real - t6_real;
            const ATOMIC t9_imag = t7_imag - t6_imag;
            
            /* t10 = - (sin(2 pi/5) t3 + sin(2 pi/10) t4) */
            const ATOMIC t10_real = -sina * t3_real - sinb * t4_real;
            const ATOMIC t10_imag = -sina * t3_imag - sinb * t4_imag;
            
            /* t11 = -(sin(2 pi/10) t3 - sin(2 pi/5) t4) */
            const ATOMIC t11_real = -sinb * t3_real + sina * t4_real;
            const ATOMIC t11_imag = -sinb * t3_imag + sina * t4_imag;
            
            /* x0 = z0 + t5 */
            const ATOMIC x0_real = z0_real + t5_real;
            const ATOMIC x0_imag = z0_imag + t5_imag;
            
            /* x1 = t8 + i t10 */
            const ATOMIC x1_real = t8_real - t10_imag;
            const ATOMIC x1_imag = t8_imag + t10_real;
            
            /* x2 = t9 + i t11 */
            const ATOMIC x2_real = t9_real - t11_imag;
            const ATOMIC x2_imag = t9_imag + t11_real;
            
            /* x3 = t9 - i t11 */
            const ATOMIC x3_real = t9_real + t11_imag;
            const ATOMIC x3_imag = t9_imag - t11_real;
            
            /* x4 = t8 - i t10 */
            const ATOMIC x4_real = t8_real + t10_imag;
            const ATOMIC x4_imag = t8_imag - t10_real;
            
            const size_t to0 = k1 * product + 2 * k - 1;
            const size_t to1 = to0 + 2 * product_1;
            const size_t to2 = to1 + 2 * product_1;
            const size_t to3 = 2 * product_1 - 2 * k + k1 * product - 1;
            const size_t to4 = to3 + 2 * product_1;
            
            VECTOR(out,ostride,to0) = x0_real;
            VECTOR(out,ostride,to0 + 1) = x0_imag;
            
            VECTOR(out,ostride,to1) = x1_real;
            VECTOR(out,ostride,to1 + 1) = x1_imag;
            
            VECTOR(out,ostride,to2) = x2_real;
            VECTOR(out,ostride,to2 + 1) = x2_imag;
            
            VECTOR(out,ostride,to3) = x4_real;
            VECTOR(out,ostride,to3 + 1) = -x4_imag;
            
            VECTOR(out,ostride,to4) = x3_real;
            VECTOR(out,ostride,to4 + 1) = -x3_imag;
        }
    }
    
    if (product_1 % 2 == 1)
        return;
    
    for (k1 = 0; k1 < q; k1++)
    {
        const size_t from0 = k1 * product_1 + product_1 - 1;
        const size_t from1 = from0 + m;
        const size_t from2 = from1 + m;
        const size_t from3 = from2 + m;
        const size_t from4 = from3 + m;
        
        const ATOMIC z0_real = VECTOR(in,istride,from0);
        const ATOMIC z1_real = VECTOR(in,istride,from1);
        const ATOMIC z2_real = VECTOR(in,istride,from2);
        const ATOMIC z3_real = VECTOR(in,istride,from3);
        const ATOMIC z4_real = VECTOR(in,istride,from4);
        
        const ATOMIC t1 = z1_real - z4_real;
        const ATOMIC t2 = z1_real + z4_real;
        const ATOMIC t3 = z2_real - z3_real;
        const ATOMIC t4 = z2_real + z3_real;
        const ATOMIC t5 = t1 - t3;
        const ATOMIC t6 = z0_real + t5 / 4.0;
        const ATOMIC t7 = (sqrt (5.0) / 4.0) * (t1 + t3);
        
        const size_t to0 = k1 * product + product_1 - 1;
        const size_t to1 = to0 + 2 * product_1;
        const size_t to2 = to1 + 2 * product_1;
        
        VECTOR(out,ostride,to0) = t6 + t7;
        VECTOR(out,ostride,to0 + 1) = -sinb * t2 - sina * t4;
        
        VECTOR(out,ostride,to1) = t6 - t7;
        VECTOR(out,ostride,to1 + 1) = -sina * t2 + sinb * t4;
        
        VECTOR(out,ostride,to2) = z0_real - t5;
    }
    
    return;

}

static void fft_real_pass_n (const double in[],
                                       const size_t istride,
                                       double out[],
                                       const size_t ostride,
                                       const size_t factor,
                                       const size_t product,
                                       const size_t n,
                             const gsl_complex twiddle[]){
    size_t k, k1;
    
    const size_t m = n / factor;
    const size_t q = n / product;
    const size_t product_1 = product / factor;
    
    size_t e1, e2;
    
    const double d_theta = 2.0 * M_PI / ((double) factor);
    const ATOMIC cos_d_theta = cos (d_theta);
    const ATOMIC sin_d_theta = sin (d_theta);
    
    for (k1 = 0; k1 < q; k1++)
    {
        /* compute x = W(factor) z, for z real */
        
        ATOMIC dw_real = 1.0, dw_imag = 0.0;
        
        for (e1 = 0; e1 <= factor - e1; e1++)
        {
            ATOMIC sum_real = 0.0;
            ATOMIC sum_imag = 0.0;
            
            ATOMIC w_real = 1.0, w_imag = 0.0;
            
            if (e1 > 0)
            {
                ATOMIC tmp_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
                ATOMIC tmp_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
                dw_real = tmp_real;
                dw_imag = tmp_imag;
            }
            
            for (e2 = 0; e2 < factor; e2++)
            {
                ATOMIC z_real = VECTOR(in,istride,k1 * product_1 + e2 * m);
                
                if (e2 > 0)
                {
                    ATOMIC tmp_real = dw_real * w_real - dw_imag * w_imag;
                    ATOMIC tmp_imag = dw_real * w_imag + dw_imag * w_real;
                    w_real = tmp_real;
                    w_imag = tmp_imag;
                }
                
                sum_real += w_real * z_real;
                sum_imag += w_imag * z_real;
                
            }
            if (e1 == 0)
            {
                const size_t to0 = product * k1;
                VECTOR(out,ostride,to0) = sum_real;
            }
            else if (e1 < factor - e1)
            {
                const size_t to0 = k1 * product + 2 * e1 * product_1 - 1;
                VECTOR(out,ostride,to0) = sum_real;
                VECTOR(out,ostride,to0 + 1) = sum_imag;
            }
            else if (e1 == factor - e1)
            {
                const size_t to0 = k1 * product + 2 * e1 * product_1 - 1;
                VECTOR(out,ostride,to0) = sum_real;
            }
            
        }
    }
    
    if (product_1 == 1)
        return;
    
    for (k = 1; k < (product_1 + 1) / 2; k++)
    {
        for (k1 = 0; k1 < q; k1++)
        {
            
            ATOMIC dw_real = 1.0, dw_imag = 0.0;
            
            for (e1 = 0; e1 < factor; e1++)
            {
                ATOMIC sum_real = 0.0, sum_imag = 0.0;
                
                ATOMIC w_real = 1.0, w_imag = 0.0;
                
                if (e1 > 0)
                {
                    const ATOMIC tmp_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
                    const ATOMIC tmp_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
                    dw_real = tmp_real;
                    dw_imag = tmp_imag;
                }
                
                for (e2 = 0; e2 < factor; e2++)
                {
                    
                    int tskip = (product_1 + 1) / 2 - 1;
                    const size_t from0 = k1 * product_1 + 2 * k + e2 * m - 1;
                    ATOMIC tw_real, tw_imag;
                    ATOMIC z_real, z_imag;
                    
                    if (e2 == 0)
                    {
                        tw_real = 1.0;
                        tw_imag = 0.0;
                    }
                    else
                    {
                        const size_t t_index = (k - 1) + (e2 - 1) * tskip;
                        tw_real = GSL_REAL(twiddle[t_index]);
                        tw_imag = -GSL_IMAG(twiddle[t_index]);
                    }
                    
                    {
                        const ATOMIC f0_real = VECTOR(in,istride,from0);
                        const ATOMIC f0_imag = VECTOR(in,istride,from0 + 1);
                        
                        z_real = tw_real * f0_real - tw_imag * f0_imag;
                        z_imag = tw_real * f0_imag + tw_imag * f0_real;
                    }
                    
                    if (e2 > 0)
                    {
                        const ATOMIC tmp_real = dw_real * w_real - dw_imag * w_imag;
                        const ATOMIC tmp_imag = dw_real * w_imag + dw_imag * w_real;
                        w_real = tmp_real;
                        w_imag = tmp_imag;
                    }
                    
                    sum_real += w_real * z_real - w_imag * z_imag;
                    sum_imag += w_real * z_imag + w_imag * z_real;
                }
                
                if (e1 < factor - e1)
                {
                    const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
                    VECTOR(out,ostride,to0) = sum_real;
                    VECTOR(out,ostride,to0 + 1) = sum_imag;
                }
                else
                {
                    const size_t to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
                    VECTOR(out,ostride,to0) = sum_real;
                    VECTOR(out,ostride,to0 + 1) = -sum_imag;
                }
                
            }
        }
    }
    
    
    if (product_1 % 2 == 1)
        return;
    
    {
        double tw_arg = M_PI / ((double) factor);
        ATOMIC cos_tw_arg = cos (tw_arg);
        ATOMIC sin_tw_arg = -sin (tw_arg);
        
        for (k1 = 0; k1 < q; k1++)
        {
            ATOMIC dw_real = 1.0, dw_imag = 0.0;
            
            for (e1 = 0; e1 < factor; e1++)
            {
                ATOMIC z_real, z_imag;
                
                ATOMIC sum_real = 0.0;
                ATOMIC sum_imag = 0.0;
                
                ATOMIC w_real = 1.0, w_imag = 0.0;
                ATOMIC tw_real = 1.0, tw_imag = 0.0;
                
                if (e1 > 0)
                {
                    ATOMIC t_real = dw_real * cos_d_theta + dw_imag * sin_d_theta;
                    ATOMIC t_imag = -dw_real * sin_d_theta + dw_imag * cos_d_theta;
                    dw_real = t_real;
                    dw_imag = t_imag;
                }
                
                for (e2 = 0; e2 < factor; e2++)
                {
                    
                    if (e2 > 0)
                    {
                        ATOMIC tmp_real = tw_real * cos_tw_arg - tw_imag * sin_tw_arg;
                        ATOMIC tmp_imag = tw_real * sin_tw_arg + tw_imag * cos_tw_arg;
                        tw_real = tmp_real;
                        tw_imag = tmp_imag;
                    }
                    
                    if (e2 > 0)
                    {
                        ATOMIC tmp_real = dw_real * w_real - dw_imag * w_imag;
                        ATOMIC tmp_imag = dw_real * w_imag + dw_imag * w_real;
                        w_real = tmp_real;
                        w_imag = tmp_imag;
                    }
                    
                    
                    {
                        const size_t from0 = k1 * product_1 + 2 * k + e2 * m - 1;
                        const ATOMIC f0_real = VECTOR(in,istride,from0);
                        z_real = tw_real * f0_real;
                        z_imag = tw_imag * f0_real;
                    }
                    
                    sum_real += w_real * z_real - w_imag * z_imag;
                    sum_imag += w_real * z_imag + w_imag * z_real;
                }
                
                if (e1 + 1 < factor - e1)
                {
                    const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
                    VECTOR(out,ostride,to0) = sum_real;
                    VECTOR(out,ostride,to0 + 1) = sum_imag;
                }
                else if (e1 + 1 == factor - e1)
                {
                    const size_t to0 = k1 * product - 1 + 2 * e1 * product_1 + 2 * k;
                    VECTOR(out,ostride,to0) = sum_real;
                }
                else
                {
                    const size_t to0 = k1 * product - 1 + 2 * (factor - e1) * product_1 - 2 * k;
                    VECTOR(out,ostride,to0) = sum_real;
                    VECTOR(out,ostride,to0 + 1) = -sum_imag;
                }
                
            }
        }
    }
    return;

}
