#ifndef MULSIG_PARAMS
#define MULSIG_PARAMS

// #define DEBUG_MULSIG

//====== ADJUSTABLE PARAMS ======
// Possible settings: 1, 2, 3, 4
#define PARAM_SET 1

//====== BE CAREFUL TO CHANGE THE FOLLOWING PARAMS  =======
#define MOD_38 274810798081L
#define MOD_39 549753716737L
#define MOD_40 1095216660481L

#define SIGMA_Y_30 1522429331L
#define SIGMA_Y_31 2309154838L
#define SIGMA_Y_32 3493703138L

#define SIGMA_SK 4 // Discrete gaussian stdev for sk
#define DEG 256    // Poynomial degree + 1
#define ROWS 6     // Rows of matrix A
#define COLS 4     // Columns
#define M_REJ 1.0  // Rejection Sampling constant factor -- M
#define TAU 23     // L1 bound of the Challenge

#if (PARAM_SET == 1)
#define MODULUS_BITS 38
#define FIXED_MODULUS (MOD_38)
#define L_COMMS 4            // # of Commitments: NEEDS TO BE <= 10
#define N_USERS 2            // # of Users
#define DELTA 31             // 2^DELTA \approx SIGMA_Y
#define SIGMA_Y (SIGMA_Y_30) // Discrete gaussian stdev for y
#elif (PARAM_SET == 2)
#define MODULUS_BITS 39
#define FIXED_MODULUS (MOD_39)
#define L_COMMS 3            // # of Commitments: NEEDS TO BE <= 10
#define N_USERS 3            // # of Users
#define DELTA 31             // 2^DELTA \approx SIGMA_Y
#define SIGMA_Y (SIGMA_Y_31) // Discrete gaussian stdev for y
#elif (PARAM_SET == 3)
#define MODULUS_BITS 40
#define FIXED_MODULUS (MOD_40)
#define L_COMMS 2            // # of Commitments: NEEDS TO BE <= 10
#define N_USERS 5            // # of Users
#define DELTA 32             // 2^DELTA \approx SIGMA_Y
#define SIGMA_Y (SIGMA_Y_32) // Discrete gaussian stdev for y
#elif (PARAM_SET == 4)
#define MODULUS_BITS 40
#define FIXED_MODULUS (MOD_40)
#define L_COMMS 2            // # of Commitments: NEEDS TO BE <= 10
#define N_USERS 7            // # of Users
#define DELTA 32             // 2^DELTA \approx SIGMA_Y
#define SIGMA_Y (SIGMA_Y_32) // Discrete gaussian stdev for y
#endif

static constexpr long MODULUS_BYTES = (MODULUS_BITS & 0x07) == 0 ? MODULUS_BITS >> 3 : (MODULUS_BITS >> 3) + 1;
#ifdef FIXED_MODULUS
#define MODULUS_FOR_INIT FIXED_MODULUS
#else
#define MODULUS_FOR_INIT (NTL::GenPrime_ZZ(MODULUS_BITS))
#endif

#define SEEDBYTES 32

static constexpr int POLYBYTES = MODULUS_BYTES * DEG;
static constexpr int OTHERS = N_USERS - 1;

static constexpr int SUM_OF_DIMS = ROWS + COLS;
static constexpr int L_MULT_ROWS = L_COMMS * ROWS;
static constexpr int L_MULT_COLS = L_COMMS * COLS;
static constexpr int L_MULT_SUM_OF_DIMS = L_COMMS * SUM_OF_DIMS;

static constexpr int TEMP_BITS = DEG * SUM_OF_DIMS * (DELTA + 13);
static constexpr int PKBYTES = DEG * ROWS * MODULUS_BITS / 8;
#endif
