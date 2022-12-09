#ifndef POLY_UTILS
#define POLY_UTILS

#include <string>
#include <assert.h>
#include <NTL/RR.h>
#include <bitset>

#include "fips202.h"
#include "sample.hpp"
#include "params.h"

inline std::string InverseIntIndex(size_t x)
{
#ifdef DEBUG_INTINDEX
    std::cout << x << " , ";
#endif
    std::string s;
    uint8_t m;
    assert(L_COMMS <= 10);
    do
    {
        m = x % L_COMMS;
        s.push_back('0' + m);
    } while (x /= L_COMMS);
    while (s.size() < N_USERS)
        s.push_back('0');
    std::reverse(s.begin(), s.end());
#ifdef DEBUG_INTINDEX
    std::cout << s << std::endl;
#endif
    return s;
}

inline int IntIndex(const int indices[])
{
    int x = 0;
    for (int i = N_USERS - 1; i >= 0; --i)
    {
        x += indices[i] * pow(L_COMMS, N_USERS - 1 - i);
    }
    return x;
}

inline void conv2Signed(NTL::ZZ &dest, const NTL::ZZ_p &src)
{
    NTL::ZZ mod = NTL::ZZ_p::modulus();
    NTL::ZZ half_mod = mod >> 1;
    NTL::ZZ temp = NTL::rep(src);
    dest = (temp < half_mod ? temp : temp - mod);
}

inline void conv2Signed(NTL::ZZX &dest, const NTL::ZZ_pX &src)
{
    assert((DEG - 1) >= NTL::deg(src));

    dest.SetLength(DEG);
    NTL::ZZ mod = NTL::ZZ_p::modulus();
    NTL::ZZ half_mod = mod >> 1;
    NTL::ZZ temp;
    for (auto i = 0; i < DEG; ++i)
    {
        temp = NTL::rep(src[i]);
        dest[i] = (temp < half_mod ? temp : temp - mod);
    }
    dest.normalize();
}

inline void conv2Signed(NTL::vec_ZZX &dest, const NTL::vec_ZZ_pX &src)
{
    dest.SetLength(src.length());
    for (auto i = 0; i < dest.length(); ++i)
    {
        conv2Signed(dest[i], src[i]);
    }
}

inline void conv2Signed(NTL::vec_ZZ &dest, const NTL::vec_ZZ_pX &src)
{

    dest.SetLength(src.length() * DEG);
    for (auto i = 0; i < src.length(); ++i)
    {
        assert((DEG - 1) >= NTL::deg(src[i]));
        for (auto j = 0; j < DEG; ++j)
        {
            conv2Signed(dest[i * DEG + j], src[i][j]);
        }
    }
}

inline void vecInnerProd(NTL::ZZ_pX &dest, const NTL::ZZ_pXModulus &modp,
                         const NTL::vec_ZZ_pX &src1, const NTL::vec_ZZ_pX &src2)
{
    assert(src1.length() == src2.length());

    NTL::ZZ_pX temp;
    auto vec_len = src1.length();
    dest = 0;
    for (auto i = 0; i < vec_len; ++i)
    {
        temp = NTL::MulMod(src1[i], src2[i], modp);
        dest += temp;
    }
}

inline void poly_pack(uint8_t packed[], const NTL::ZZ_pX &poly)
{
    assert((DEG - 1) >= NTL::deg(poly));

    NTL::ZZ tmp;
    for (auto i = 0; i < DEG; ++i)
    {
        tmp = rep(poly[i]);
        NTL::BytesFromZZ(packed + i * MODULUS_BYTES, tmp, MODULUS_BYTES);
    }
}

inline void polyvec_pack(uint8_t packed[],
                         const NTL::vec_ZZ_pX &polyvec)
{
    auto len = polyvec.length();

    for (auto i = 0; i < len; ++i)
    {
        poly_pack(packed + i * POLYBYTES, polyvec[i]);
    }
}

inline void polyvec_pack(uint8_t packed[],
                         const NTL::vec_ZZ_pX &polyvec, int from, int len)
{
    for (auto i = 0; i < len; ++i)
    {
        poly_pack(packed + i * POLYBYTES, polyvec[from + i]);
    }
}

inline void polyvec_pack(uint8_t packed[],
                         const NTL::Vec<NTL::vec_ZZ_pX> &polyvec)
{
    auto len = polyvec.length();
    auto veclen = polyvec[0].length();
    for (auto i = 0; i < len; ++i)
    {
        assert(polyvec[i].length() == veclen);
        polyvec_pack(packed + i * POLYBYTES * veclen, polyvec[i]);
    }
}

inline void poly_unpack(NTL::ZZ_pX &poly, const uint8_t packed[])
{
    NTL::ZZ temp;

    poly.SetLength(DEG);
    for (auto i = 0; i < DEG; ++i)
    {
        temp = NTL::ZZFromBytes(packed, DEG);
        poly[i] = NTL::conv<NTL::ZZ_p>(temp);
    }
    poly.normalize();
}

inline void polyvec_unpack(NTL::vec_ZZ_pX &polyvec,
                           const uint8_t packed[], int veclen)
{
    polyvec.SetLength(veclen);

    for (auto i = 0; i < veclen; ++i)
    {
        poly_unpack(polyvec[i], packed + i * POLYBYTES);
    }
}

// Samples polynomial with TAU nonzero
// coefficients in {-1,1} using the output stream of
// SHAKE256(seed).
inline void poly_challenge(NTL::ZZ_pX &c,
                           const uint8_t seed[SEEDBYTES])
{
    assert(DEG >= TAU);

    unsigned int i, b, pos;
    uint64_t signs;
    uint8_t buf[SHAKE256_RATE];
    keccak_state state;

    c.SetLength(DEG);

    shake256_init(&state);
    shake256_absorb(&state, seed, SEEDBYTES);
    shake256_finalize(&state);
    shake256_squeezeblocks(buf, 1, &state);

    signs = 0;
    for (i = 0; i < 8; ++i)
        signs |= (uint64_t)buf[i] << 8 * i;
    pos = 8;

    for (i = DEG - TAU; i < DEG; ++i)
    {
        do
        {
            if (pos >= SHAKE256_RATE)
            {
                shake256_squeezeblocks(buf, 1, &state);
                pos = 0;
            }

            b = buf[pos++];
        } while (b > i);

        c[i] = c[b];
        c[b] = 1 - 2 * (signs & 1);
        signs >>= 1;
    }

    c.normalize();
}

inline bool rejSampOnce(const NTL::vec_ZZ &sample,
                        const NTL::vec_ZZ &center)
{
    bool res = true;
    NTL::ZZ centerSq, innerProd, temp;
    NTL::RR exponent, prob;
    NTL::RR sigmaSq{double(SIGMA_Y)};

    assert(sample.length() == center.length());

    NTL::InnerProduct(innerProd, sample, center);
    NTL::InnerProduct(centerSq, center, center);

    NTL::mul(temp, innerProd, -2);
    NTL::add(temp, temp, centerSq);

    NTL::sqr(sigmaSq, sigmaSq);
    NTL::mul(sigmaSq, sigmaSq, 2);

    exponent = NTL::conv<NTL::RR>(temp);
    NTL::div(exponent, exponent, sigmaSq);

    NTL::exp(prob, exponent);
    NTL::div(prob, prob, double(M_REJ));

    NTL::RR thresh = NTL::random_RR();

#ifdef DEBUG_MULSIG
    std::cout << "prob: " << prob;
    std::cout << ", thresh: " << thresh << std::endl;
#endif
    if (prob >= thresh) // accept
        res = true;
    else
        res = false;

    return res;
}

inline int rejSamp(const NTL::vec_ZZ_pX &y,
                   const NTL::vec_ZZ_pX &z)
{
    int acc_idx = -1;
    NTL::vec_ZZ_pX x;
    NTL::vec_ZZ z_signed, x_signed;
    x.SetLength(SUM_OF_DIMS);

    conv2Signed(z_signed, z);

    for (auto l = 0; l < L_COMMS; ++l)
    {
        for (auto i = 0; i < SUM_OF_DIMS; ++i)
        {
            x[i] = y[i + l * SUM_OF_DIMS] + z[i];
        }
        conv2Signed(x_signed, x);
        if (rejSampOnce(x_signed, z_signed))
        {
            acc_idx = l;
            break;
        }
    }

    return acc_idx;
}

inline int encode(std::bitset<TEMP_BITS> &code, const NTL::ZZ_pX &z)
{
    assert((DEG - 1) >= NTL::deg(z));
    NTL::ZZ signed_z, z1;
    long delt = 1L << DELTA;
    int idx = 0;
    for (auto d = 0; d < DEG; ++d)
    {
        int k = 0, num_zeros = 0;
        long z0 = 0;
        unsigned long uz0 = 0;
        conv2Signed(signed_z, z[d]);
        z0 = NTL::DivRem(z1, signed_z, delt);
        if (z0 >= (delt >> 1))
        {
            uz0 = delt - z0;
            code[idx] = 1;
        }
        else
        {
            uz0 = z0;
            code[idx] = 0;
        }
        idx++;
#ifdef DEBUG_MULSIG_ENCODE
        std::cout << uz0 << ", ";
#endif
        for (auto i = 0; i < (DELTA - 1); ++i)
        {
            code[idx] = (uz0 & 0x01);
            uz0 >>= 1;
            idx++;
        }
        k = NTL::conv<int>(z1);
        switch (k)
        {
        case 0:
            code[idx] = 0;
            code[idx + 1] = 0;
            break;
        case 1:
            code[idx] = 0;
            code[idx + 1] = 1;
            break;
        case -1:
            code[idx] = 1;
            code[idx + 1] = 0;
            break;
        default:
            code[idx] = 1;
            code[idx + 1] = 1;
            break;
        }
        idx += 2;
        if (code[idx - 1] == 1 && code[idx - 2] == 1)
        {
            num_zeros = (k > 0) ? 2 * k - 4 : 2 * (-k) - 3;
            for (auto i = 0; i < num_zeros; ++i)
            {
                code[idx] = 0;
                idx++;
            }
            code[idx] = 1;
            idx++;
        }
    }
    return idx - 1; // number of bits used
}
inline int encode_polyvec(NTL::Vec<std::bitset<TEMP_BITS>> &vec_code,
                          const NTL::vec_ZZ_pX &vec_z)
{
    auto veclen = vec_z.length();
    vec_code.SetLength(veclen);
    int total_bits = 0;
    for (auto i = 0; i < veclen; ++i)
    {
        total_bits += encode(vec_code[i], vec_z[i]);
    }
    return total_bits;
}

inline void decode(NTL::ZZ_pX &z, const std::bitset<TEMP_BITS> &code)
{
    z.SetLength(DEG);
    NTL::ZZ z1, signed_z;
    long delt = 1L << DELTA;
    int idx = 0;
    for (auto d = 0; d < DEG; ++d)
    {
        int num_zeros = 0;
        long z0 = 0;
        unsigned long uz0 = 0;
        int sign = code[idx];
        idx++;
        for (auto i = 0; i < (DELTA - 1); ++i)
        {
            uz0 |= (code[idx] << i);
            idx++;
        }
        z0 = (sign == 0) ? uz0 : delt - uz0;
#ifdef DEBUG_MULSIG_DECODE
        std::cout << uz0 << ", ";
#endif
        if (code[idx] == 0 && code[idx + 1] == 0)
        {
            z1 = 0;
            idx += 2;
        }
        else if (code[idx] == 0 && code[idx + 1] == 1)
        {
            z1 = 1;
            idx += 2;
        }
        else if (code[idx] == 1 && code[idx + 1] == 0)
        {
            z1 = -1;
            idx += 2;
        }
        else
        { // 0x11
            idx += 2;
            while (code[idx] == 0)
            {
                num_zeros++;
                idx++;
            }
            if (code[idx] != 1)
                throw std::runtime_error("Error: Huffman decoding failed!");
            z1 = (num_zeros & 0x01) == 0 ? (num_zeros + 4) / 2 : -((num_zeros + 3) / 2);
            idx++;
        }
        signed_z = z1 * delt + z0;
        z[d] = NTL::conv<NTL::ZZ_p>(signed_z);
    }
}

inline void decode_polyvec(NTL::vec_ZZ_pX &vec_z,
                           const NTL::Vec<std::bitset<TEMP_BITS>> &vec_code)
{
    auto veclen = vec_code.length();
    vec_z.SetLength(veclen);
    for (auto i = 0; i < veclen; ++i)
    {
        decode(vec_z[i], vec_code[i]);
    }
}

#endif
