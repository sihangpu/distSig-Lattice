#ifndef MULSIG_H
#define MULSIG_H

#include <random>
#include <chrono>
#include "merklecpp.h"
#include "utils.hpp"

inline void init(NTL::ZZ_pXModulus &modp, NTL::Vec<NTL::vec_ZZ_pX> &matA)
{
    NTL::ZZ _mod(MODULUS_FOR_INIT);
    NTL::ZZ_p::init(_mod);
#ifdef DEBUG_MULSIG
    std::cout << "Modulus in usage: " << NTL::ZZ_p::modulus() << std::endl;
#endif
    // cyclotomic polynomial: phi_{2*DEG}(X)= X^DEG+1
    NTL::ZZ_pX _poly = NTL::ZZ_pX(NTL::INIT_MONO, DEG);
    NTL::SetCoeff(_poly, 0, 1);
    NTL::build(modp, _poly);

    matA.SetLength(ROWS);
    for (auto i = 0; i < ROWS; ++i)
    {
        matA[i].SetLength(COLS);
        for (auto j = 0; j < COLS; ++j)
        {
            matA[i][j] = NTL::random_ZZ_pX(DEG);
        }
    }
}

inline void key_pair(NTL::vec_ZZ_pX &pk, NTL::vec_ZZ_pX &sk,
                     const NTL::Vec<NTL::vec_ZZ_pX> &A,
                     const NTL::ZZ_pXModulus &modp)
{
    NTL::xdouble sigma_sk(double(SIGMA_SK));
    NTL::ZZ_pX temp;

    sk.SetLength(SUM_OF_DIMS);
    pk.SetLength(ROWS);

    for (auto i = 0; i < COLS; ++i)
    {
        sampleGaussian(sk[i], DEG, sigma_sk);
    }

    for (auto i = 0; i < ROWS; ++i)
    {
        sampleGaussian(sk[i + COLS], DEG, sigma_sk);
        pk[i] = sk[i + COLS];
        for (auto j = 0; j < COLS; ++j)
        {
            temp = NTL::MulMod(A[i][j], sk[j], modp);
            pk[i] += temp;
        }
    }

#ifdef DEBUG_MULSIG_KEYGEN
    std::cout << "PK:\n"
              << pk << std::endl
              << "SK:\n"
              << sk << std::endl;

#endif
}

inline void sign_1st(NTL::vec_ZZ_pX &y, NTL::vec_ZZ_pX &v,
                     const NTL::Vec<NTL::vec_ZZ_pX> &A,
                     const NTL::ZZ_pXModulus &modp)
{
    int step_y, step_v;
    NTL::xdouble sigma_y(double(SIGMA_Y));
    NTL::ZZ_pX temp;

    y.SetLength(L_MULT_SUM_OF_DIMS);
    v.SetLength(L_MULT_ROWS);

    for (auto l = 0; l < L_COMMS; ++l)
    {
        step_y = l * SUM_OF_DIMS;
        step_v = l * ROWS;
        for (auto i = 0; i < SUM_OF_DIMS; ++i)
        {
            sampleGaussian(y[step_y + i], DEG, sigma_y);
        }
        for (auto i = 0; i < ROWS; ++i)
        {
            v[step_v + i] = y[step_y + COLS + i];
            for (auto j = 0; j < COLS; ++j)
            {
                temp = NTL::MulMod(A[i][j], y[step_y + j], modp);
                v[step_v + i] += temp;
            }
        }
    }
#ifdef DEBUG_MULSIG_SIGN_1
    std::cout << "y:\n"
              << y << std::endl
              << "v:\n"
              << v << std::endl;
#endif
}

inline int sign_2nd(NTL::vec_ZZ_pX &z, int &rejSamp_id, merkle::Tree &mt, int user_id,
                    const NTL::Vec<NTL::vec_ZZ_pX> &A, const NTL::ZZ_pXModulus &modp,
                    const NTL::Vec<NTL::vec_ZZ_pX> &pk_list, const NTL::vec_ZZ_pX &sk,
                    const NTL::vec_ZZ_pX &y, const NTL::Vec<NTL::vec_ZZ_pX> &v_list,
                    const uint8_t msg[], size_t mlen)
{
    // @TODO: check integrity of v_list
    // then do the following
    size_t total_leaves = pow(L_COMMS, N_USERS);
    NTL::vec_ZZ_pX w;
    w.SetLength(total_leaves * ROWS);
    z.SetLength(SUM_OF_DIMS);

    int idx;
    // Get w from all users' v
    for (size_t i = 0; i < total_leaves; ++i)
    {
        std::string tempstr = InverseIntIndex(i);
        if (tempstr.size() != N_USERS)
            throw std::runtime_error("ERR:sign_2nd: Errors in InvIntIndex!");

        for (size_t j = 0; j < N_USERS; ++j)
        {
            idx = std::stoi(tempstr.substr(j, 1));
            for (auto k = 0; k < ROWS; ++k)
            {
                w[i * ROWS + k] += v_list[j][idx * ROWS + k];
            }
        }
    }
    // Build the merkle tree
    std::vector<merkle::Hash> hashes(total_leaves);
    uint8_t binary_w[ROWS * POLYBYTES];
    uint8_t hash[SEEDBYTES];
    for (auto i = 0; i < total_leaves; ++i)
    {
        polyvec_pack(binary_w, w, i * ROWS, ROWS);
        shake256(hash, SEEDBYTES, binary_w, sizeof(binary_w));
        hashes[i] = hash;
    }
    for (auto h : hashes)
        mt.insert(h);

    auto root = mt.root();

    // Get the challenge
    uint8_t seed_for_c[SEEDBYTES];
    uint8_t binary_pk_rt_m[POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES + mlen];
    NTL::ZZ_pX c;
    polyvec_pack(binary_pk_rt_m, pk_list[user_id]);
    polyvec_pack(binary_pk_rt_m + POLYBYTES * ROWS, pk_list);
    memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS, root.bytes, SEEDBYTES);
    memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES, msg, mlen);
    shake256(seed_for_c, SEEDBYTES, binary_pk_rt_m, sizeof(binary_pk_rt_m));
    poly_challenge(c, seed_for_c);
#ifdef DEBUG_MULSIG_PACK
    std::cout << pk_list[user_id] << std::endl;
    for (auto byte : binary_pk_rt_m)
        printf("%d, ", byte);
    std::cout << c << std::endl;
#endif
    for (auto i = 0; i < SUM_OF_DIMS; ++i)
    {
        NTL::MulMod(z[i], sk[i], c, modp);
    }

    // Rejection sampling
    rejSamp_id = rejSamp(y, z);

    if ((rejSamp_id < 0) || (rejSamp_id >= L_COMMS))
    {
#ifdef DEBUG_MULSIG
        printf("RejSampling for User #%d: Failed.\n", user_id);
#endif
        return 0xFF;
    }
    for (auto i = 0; i < SUM_OF_DIMS; ++i)
    {
        z[i] += y[rejSamp_id * SUM_OF_DIMS + i];
    }

#ifdef DEBUG_MULSIG
    printf("RejSampling for User #%d: Passed.\n", user_id);
#endif
    return 0;
}

inline bool verify_id(int user_id, int &rejSamp_id,
                      const NTL::Vec<NTL::vec_ZZ_pX> &A, const NTL::ZZ_pXModulus &modp,
                      const NTL::Vec<NTL::vec_ZZ_pX> &pk_list, const NTL::vec_ZZ_pX &v,
                      const NTL::vec_ZZ_pX &z, const merkle::Hash &root,
                      const uint8_t msg[], size_t mlen)
{
    // @TODO: check ||z|| > B_z
    // then do the following
    // hash and get the challenge
    uint8_t seed_for_c[SEEDBYTES];
    uint8_t binary_pk_rt_m[POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES + mlen];
    NTL::ZZ_pX c;
    polyvec_pack(binary_pk_rt_m, pk_list[user_id]);
    polyvec_pack(binary_pk_rt_m + POLYBYTES * ROWS, pk_list);
    memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS, root.bytes, SEEDBYTES);
    memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES, msg, mlen);
    shake256(seed_for_c, SEEDBYTES, binary_pk_rt_m, sizeof(binary_pk_rt_m));
    poly_challenge(c, seed_for_c);
#ifdef DEBUG_MULSIG_PACK
    std::cout << pk_list[user_id] << std::endl;
    for (auto byte : binary_pk_rt_m)
        printf("%d, ", byte);
    std::cout << c << std::endl;
#endif
    NTL::vec_ZZ_pX w;
    NTL::ZZ_pX temp;
    uint8_t equal = 0x01;
    w.SetLength(ROWS);

    for (auto i = 0; i < ROWS; ++i)
    {
        NTL::MulMod(temp, pk_list[user_id][i], c, modp);
        w[i] = z[i + COLS] - temp;
        for (auto j = 0; j < COLS; ++j)
        {
            NTL::MulMod(temp, A[i][j], z[j], modp);
            w[i] += temp;
        }
    }

    for (auto i = 0; i < L_COMMS; ++i)
    {
        for (auto j = 0; j < ROWS; ++j)
        {
            equal &= (w[j] == v[i * ROWS + j]);
        }
        if (equal == 0x01)
        {
            rejSamp_id = i;
            return true;
        }
        equal = 0x01;
    }
    return false;
}

inline void sign_3rd(merkle::Hash &sig_root, std::shared_ptr<merkle::Path> &sig_path_ptr,
                     NTL::vec_ZZ_pX &sig_z, int user_id, const NTL::Vec<NTL::vec_ZZ_pX> &A,
                     const NTL::ZZ_pXModulus &modp, const NTL::Vec<NTL::vec_ZZ_pX> &pk_list,
                     const NTL::Vec<NTL::vec_ZZ_pX> &z_list, const NTL::Vec<NTL::vec_ZZ_pX> &v_list,
                     merkle::Tree mt_inputs[], const uint8_t msg[], size_t mlen, const int rejSampIds[])
{
    sig_z.SetLength(SUM_OF_DIMS);
    // Verify others' z
    int indices[N_USERS] = {-1};
    for (int i = 0; i < N_USERS; ++i)
    {
        if (i == user_id)
            continue;
        if (!verify_id(i, indices[i], A, modp, pk_list, v_list[i], z_list[i],
                       mt_inputs[i].root(), msg, mlen))
        {
            throw std::runtime_error("ERR:sign_3rd: Verify other's z failed!");
        }
        if (indices[i] != rejSampIds[i])
        {
            throw std::runtime_error("ERR:sign_3rd: rejSampIndex does not match!");
        }
    }

    for (auto i = 0; i < N_USERS; ++i)
    {
        for (auto j = 0; j < SUM_OF_DIMS; ++j)
        {
            sig_z[j] += z_list[i][j];
        }
    }

    int path_idx = IntIndex(rejSampIds);
    sig_root = mt_inputs[user_id].root();
    mt_inputs[user_id].flush_to(path_idx);
    sig_path_ptr = mt_inputs[user_id].path(path_idx);
}

inline bool verify(const NTL::Vec<NTL::vec_ZZ_pX> &A, const NTL::ZZ_pXModulus &modp,
                   const NTL::Vec<NTL::vec_ZZ_pX> &pk_list, const uint8_t msg[], int mlen,
                   const merkle::Hash &sig_root, merkle::Path &sig_path, const NTL::vec_ZZ_pX &sig_z)
{
    NTL::vec_ZZ_pX c_list, w;
    c_list.SetLength(N_USERS);
    w.SetLength(ROWS);

    NTL::ZZ_pX temp;

    uint8_t seed_for_c[SEEDBYTES];
    uint8_t binary_pk_rt_m[POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES + mlen];
    for (auto i = 0; i < N_USERS; ++i)
    {
        polyvec_pack(binary_pk_rt_m, pk_list[i]);
        polyvec_pack(binary_pk_rt_m + POLYBYTES * ROWS, pk_list);
        memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS, sig_root.bytes, SEEDBYTES);
        memcpy(binary_pk_rt_m + POLYBYTES * (N_USERS + 1) * ROWS + SEEDBYTES, msg, mlen);
        shake256(seed_for_c, SEEDBYTES, binary_pk_rt_m, sizeof(binary_pk_rt_m));
        poly_challenge(c_list[i], seed_for_c);

        for (auto j = 0; j < ROWS; ++j)
        {
            NTL::MulMod(temp, c_list[i], pk_list[i][j], modp);
            w[j] += temp;
        }
    }
    for (auto i = 0; i < ROWS; ++i)
    {
        w[i] = sig_z[i + COLS] - w[i];
        for (auto j = 0; j < COLS; ++j)
        {
            NTL::MulMod(temp, A[i][j], sig_z[j], modp);
            w[i] += temp;
        }
    }
    uint8_t binary_w[ROWS * POLYBYTES];
    uint8_t hash[SEEDBYTES];
    polyvec_pack(binary_w, w);
    shake256(hash, SEEDBYTES, binary_w, sizeof(binary_w));
    merkle::Hash leaf_w = hash;

    if (leaf_w != sig_path.leaf())
        return false;
    if (!sig_path.verify(sig_root))
        return false;

    return true;
}

#endif
