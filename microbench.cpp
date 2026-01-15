#include <cassert>
#include <iostream>
#include <vector>
#include <chrono>

#include "poly_arith.h"
#include "ringsnark/zk_proof_systems/rinocchio/rinocchio.hpp"
#include "seal/seal.h"

using namespace std;
using namespace seal;
using namespace polytools;

#define LOG_T 30
#define SEC_PARAM 128
#define NUM_REPEATS 100

struct ParamsConfig {
    size_t N;
    vector<int> coeff_bits;
};

int main() {
    vector<ParamsConfig> configs = {
        {16384, {50, 50, 50, 50, 50, 50}},
        {4096,  {50, 50}},
        {8192,  {50, 50}},
        {16384, {50, 50}}
    };

    for (auto& cfg : configs) {
        cout << "===========================" << endl;
        cout << "[CONFIG] N=" << cfg.N << ", CoeffBits={";
        for (size_t i = 0; i < cfg.coeff_bits.size(); i++) {
            if (i != 0) cout << ",";
            cout << cfg.coeff_bits[i];
        }
        cout << "}" << endl;

        // --------------------------
        // Setup SEAL
        // --------------------------
        EncryptionParameters parms(scheme_type::bgv);
        parms.set_poly_modulus_degree(cfg.N);
        parms.set_coeff_modulus(CoeffModulus::Create(cfg.N, cfg.coeff_bits));
        parms.set_plain_modulus(PlainModulus::Batching(cfg.N, LOG_T));
        SEALContext context(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);

        Encryptor encryptor(context, public_key);
        Evaluator evaluator(context);
        Decryptor decryptor(context, secret_key);
        BatchEncoder batch_encoder(context);

        size_t slot_count = batch_encoder.slot_count();
        size_t row_size = slot_count / 2;

        // --------------------------
        // Encrypt zeros
        // --------------------------
        auto start = chrono::high_resolution_clock::now();
        vector<Ciphertext> zeros_encrypted(SEC_PARAM);
        vector<uint64_t> zeros(slot_count, 0ULL);
        Plaintext zero_plain;
        batch_encoder.encode(zeros, zero_plain);
        for (int i = 0; i < SEC_PARAM; i++) {
            encryptor.encrypt(zero_plain, zeros_encrypted[i]);
        }
        auto end = chrono::high_resolution_clock::now();
        cout << "[TIME][CLIENT] One-time setup: "
             << chrono::duration_cast<chrono::microseconds>(end - start).count()
             << " us" << endl;

        // --------------------------
        // Prepare sample vector
        // --------------------------
        vector<uint64_t> pod_matrix(slot_count);
        for (size_t i = 0; i < row_size; i++) {
            pod_matrix[i] = i;
            pod_matrix[row_size + i] = 3 * i;
        }

        Plaintext x_plain;
        batch_encoder.encode(pod_matrix, x_plain);
        Ciphertext x_encrypted;
        encryptor.encrypt(x_plain, x_encrypted);

        SealPoly x1(context, x_encrypted, 0);
        SealPoly x2(context, x_encrypted, 1);
        auto tables = context.get_context_data(x1.get_parms_id())->small_ntt_tables();

        if (x1.is_ntt_form()) {
            x1.intt_inplace(tables);
        }
        x1.ntt_inplace(tables);
        if (!x2.is_ntt_form()) x2.ntt_inplace(tables);

        // --------------------------
        // Microbenchmark timers
        // --------------------------
        size_t time_NTT = 0;
        size_t time_AxR = 0;
        size_t time_RpR = 0;
        size_t time_RxR = 0;
        size_t time_enc = 0;
        size_t time_dec = 0;

        // Separate context for enc/dec
        EncryptionParameters parms2(scheme_type::bgv);
        parms2.set_poly_modulus_degree(cfg.N);
        parms2.set_coeff_modulus(CoeffModulus::Create(cfg.N, cfg.coeff_bits));
        parms2.set_plain_modulus(PlainModulus::Batching(cfg.N, 60));
        SEALContext context2(parms2);
        KeyGenerator keygen2(context2);
        SecretKey secret_key2 = keygen2.secret_key();
        PublicKey public_key2;
        keygen2.create_public_key(public_key2);
        Encryptor encryptor2(context2, public_key2);
        Decryptor decryptor2(context2, secret_key2);
        BatchEncoder batch_encoder2(context2);

        for (int i = 0; i < NUM_REPEATS; i++) {
            // NTT
            start = chrono::high_resolution_clock::now();
            x1.ntt_inplace(tables);
            end = chrono::high_resolution_clock::now();
            time_NTT += chrono::duration_cast<chrono::microseconds>(end - start).count();
            x1.intt_inplace(tables);

            // AxR
            start = chrono::high_resolution_clock::now();
            x1.multiply_scalar_inplace(3);
            end = chrono::high_resolution_clock::now();
            time_AxR += chrono::duration_cast<chrono::microseconds>(end - start).count();

            // R + R
            start = chrono::high_resolution_clock::now();
            x1.add_inplace(x2);
            end = chrono::high_resolution_clock::now();
            time_RpR += chrono::duration_cast<chrono::microseconds>(end - start).count();

            // R x R
            start = chrono::high_resolution_clock::now();
            x1.multiply_inplace(x2);
            end = chrono::high_resolution_clock::now();
            time_RxR += chrono::duration_cast<chrono::microseconds>(end - start).count();

            // 1 Enc
            Plaintext ptxt;
            Ciphertext ctxt;
            start = chrono::high_resolution_clock::now();
            batch_encoder2.encode(pod_matrix, ptxt);
            encryptor2.encrypt(ptxt, ctxt);
            end = chrono::high_resolution_clock::now();
            time_enc += chrono::duration_cast<chrono::microseconds>(end - start).count();

            // 1 Dec
            start = chrono::high_resolution_clock::now();
            decryptor2.decrypt(ctxt, ptxt);
            batch_encoder2.decode(ptxt, pod_matrix);
            end = chrono::high_resolution_clock::now();
            time_dec += chrono::duration_cast<chrono::microseconds>(end - start).count();
        }

        // --------------------------
        // Print microbenchmark results
        // --------------------------
        cout << "[TIME] NTT:   " << float(time_NTT)/NUM_REPEATS << " us" << endl;
        cout << "[TIME] A x R: " << float(time_AxR)/NUM_REPEATS << " us" << endl;
        cout << "[TIME] R + R: " << float(time_RpR)/NUM_REPEATS << " us" << endl;
        cout << "[TIME] R x R: " << float(time_RxR)/NUM_REPEATS << " us" << endl;
        cout << "[TIME] 1 Enc: " << float(time_enc)/NUM_REPEATS << " us" << endl;
        cout << "[TIME] 1 Dec: " << float(time_dec)/NUM_REPEATS << " us" << endl;

        unsigned long long size = 9 * x_encrypted.size() *
                                  parms.coeff_modulus().size() *
                                  parms.poly_modulus_degree() * 8;
        cout << "[SPACE] Proof size: " << size << " B" << endl;
        cout << endl;
    }
}
