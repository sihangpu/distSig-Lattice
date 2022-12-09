#include <random>
#include <climits>
#include <algorithm>
#include <functional>
#include <chrono>

#include "../params.h"
#include "../mulsig.hpp"
#include "../fips202.h"

using namespace NTL;
using random_bytes_engine = std::independent_bits_engine<
    std::default_random_engine, CHAR_BIT, unsigned char>;
random_bytes_engine rbe;

Vec<std::bitset<TEMP_BITS>> sig_z_bin;
int TESTIMES = 100;

static int run()
{
    double init_st, init_ed, keypair_ed,
        sign1_ed, sign2_ed, sign3_ed, verify_ed;
    double init_t = 0.0, keypair_t = 0.0, sign_t = 0.0, verify_t = 0.0;
    int bits, pathbytesNum;
    for (auto testid = 0; testid < TESTIMES; ++testid)
    {
        ZZ_pXModulus modp;
        Vec<vec_ZZ_pX> matA;
        Vec<vec_ZZ_pX> pk_list, sk_list, y_list, v_list, z_list;
        NTL::vec_ZZ_pX sig_z;

        merkle::Tree merkleTrees[N_USERS];
        merkle::Hash sig_root;
        std::shared_ptr<merkle::Path> sig_path_ptr;
        vec_ZZ_pX z_rec;
        ZZ_pX res;
        sig_z_bin.kill();

        int rejSampIds[N_USERS];
        pk_list.SetLength(N_USERS);
        sk_list.SetLength(N_USERS);
        y_list.SetLength(N_USERS);
        v_list.SetLength(N_USERS);
        z_list.SetLength(N_USERS);

        uint8_t msg[56];
        std::generate(msg, msg + sizeof(msg), std::ref(rbe));

        init_st = GetTime();
        init(modp, matA);
        init_ed = GetTime();

        for (auto i = 0; i < N_USERS; ++i)
            key_pair(pk_list[i], sk_list[i], matA, modp);
        keypair_ed = GetTime();

        for (auto i = 0; i < N_USERS; ++i)
            sign_1st(y_list[i], v_list[i], matA, modp);
        sign1_ed = GetTime();

        int restart = 0x00;
        for (auto i = 0; i < N_USERS; ++i)
            restart |= sign_2nd(z_list[i], rejSampIds[i], merkleTrees[i], i, matA, modp,
                                pk_list, sk_list[i], y_list[i], v_list, msg, sizeof(msg));
        sign2_ed = GetTime();

#ifdef DEBUG_MULSIG_MERKLE_TREE
        for (auto mtree : merkleTrees)
        {
            std::cout << mtree.to_string(3) << std::endl;
        }
#endif
        if (restart == 0x00)
        {
            sign_3rd(sig_root, sig_path_ptr, sig_z, 0,
                     matA, modp, pk_list, z_list, v_list, merkleTrees,
                     msg, sizeof(msg), rejSampIds);
            sign3_ed = GetTime();
            if (verify(matA, modp, pk_list, msg, sizeof(msg), sig_root, *sig_path_ptr, sig_z))
#ifdef DEBUG_MULSIG
                printf(">>> Verified <<<\n");
#else
                ;
#endif
            else
                printf("<<< Verification Failed >>>\n");
            verify_ed = GetTime();
        }
        else
        {
            printf("... Restart ...\n");
        }

        bits = encode_polyvec(sig_z_bin, sig_z);
        decode_polyvec(z_rec, sig_z_bin);
        pathbytesNum = ceil(log2(pow(L_COMMS, N_USERS)));

        for (auto i = 0; i < sig_z.length(); ++i)
        {
            res += sig_z[i] - z_rec[i];
        }
        if (res != 0)
        {
            printf("Err: decode/encode failed!\n");
            return 1;
        }

        auto sign1sec = (sign1_ed - keypair_ed) / N_USERS;
        auto sign2sec = (sign2_ed - sign1_ed) / N_USERS;
        auto sign3sec = sign3_ed - sign2_ed;
        init_t += (init_ed - init_st);
        keypair_t += ((keypair_ed - init_ed) / N_USERS);
        sign_t += (sign1sec + sign2sec + sign3sec);
        verify_t += (verify_ed - sign3_ed);
    }
    printf("==========================================================\n");
    printf("Finished: Run %d times successfully\n", TESTIMES);
    printf("public key size %f kilobytes\n", PKBYTES / 1024.0);
    printf("signature size: %f kilobytes\n", (bits / 8 + SEEDBYTES * pathbytesNum) / 1024.0);
    printf("init takes %f seconds\n", init_t / TESTIMES);
    printf("keygen takes %f seconds per user\n", keypair_t / TESTIMES);
    printf("signing takes %f seconds per user\n", sign_t / TESTIMES);
    printf("verifying takes %f seconds per user\n", verify_t / TESTIMES);

    return 0;
}

int main(int argc, const char *argv[])
{

    return run();
}