import math

# --------------------------------------------------
# Circuit size derivation
# --------------------------------------------------

def compute_circuit_sizes(M, N, m):
    D = N

    G  = M * (4*N + 3*D) / m
    W  = 2*D + M * (4*N + 3*D + 1) / m

    Gp = 2*(2*D*m + M) + 44*m
    Wp = 16*D*m + 10*M + 4*D + 57*m

    logG  = math.log2(G)
    logGp = math.log2(Gp)

    return G, W, Gp, Wp, logG, logGp


# --------------------------------------------------
# Runtime functions
# --------------------------------------------------

def runtime_sub_prover(G, W,
    pi_1_ntt, pi_1_mul, pi_1_add,
    pi_2_ntt, pi_2_mul, pi_2_add,
    h_fft, base_fft, base_mul, base_add,
    field_mul=36.6, field_add=6.6, num_core=8):

    term1 = 2 * W * (pi_1_ntt/2 + pi_1_mul + pi_1_add) / num_core

    term2 = (
        h_fft * ((base_fft + base_mul) / field_mul + base_add / field_add)
        / 2 / 4
        + 2 * G * (pi_2_ntt + pi_2_mul + pi_2_add) / num_core
    )

    return term1 + term2


def runtime_master_prover(Gp, Wp,
    pi_2_ntt, pi_2_mul, pi_2_add,
    h_fft_ac, base_fft_ac, base_mul_ac, base_add_ac,
    field_mul=36.6, field_add=6.6, num_core=8):

    term1 = 2 * Wp * (pi_2_ntt/2 + pi_2_mul + pi_2_add) / num_core

    term2 = (
        h_fft_ac * ((base_fft_ac + base_mul_ac) / field_mul
                    + base_add_ac / field_add)
        / 2 / 2
        + 2 * Gp * (pi_2_ntt + pi_2_mul + pi_2_add) / num_core
    )

    return term1 + term2


def runtime_verifier(N, M,
    pi_2_ntt, pi_2_mul, pi_2_add,
    pi_2_decrypt, pi_2_pl_mul,
    num_core=8):

    D = N
    return (
        2*(2*D + M) * (pi_2_ntt + pi_2_mul + pi_2_add)
        + 11*pi_2_decrypt
        + 9*pi_2_pl_mul
        + 2*pi_2_mul
    ) / num_core


# --------------------------------------------------
# Benchmarked constants (placeholders)
# --------------------------------------------------

# ring dim: 4096, q : 109 bit
pi_1_ntt = 0.00003884
pi_1_mul = 0.00001139
pi_1_add = 0.00000301

# ring dim: 16384, Q : 300 bit
pi_2_ntt = 0.00094957
pi_2_mul = 0.00017846
pi_2_add = 0.00005294
pi_2_decrypt = 0.00158994
pi_2_pl_mul = 0.00007477

# FFT parameters stored in dictionaries
fft_params_M = {
    # 1024:  {'h_fft': 0.137299, 'h_fft_ac': 0.163174, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975, 'pi_1_ntt': 0.00003884, 'pi_1_mul': 0.00001139, 'pi_1_add': 0.00000301},
    4096:  {'h_fft': 0.301866, 'h_fft_ac': 0.163766, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975, 'pi_1_ntt': 0.00003884, 'pi_1_mul': 0.00001139, 'pi_1_add': 0.00000301},
    8192:  {'h_fft': 0.623237, 'h_fft_ac': 0.348989, 'base_fft': 1293000, 'base_mul': 300156, 'base_add': 53950, 'pi_1_ntt': 0.00008578, 'pi_1_mul': 0.00002323, 'pi_1_add': 0.00000726},
    16384: {'h_fft': 1.267, 'h_fft_ac': 0.717230, 'base_fft': 2189000, 'base_mul': 600313, 'base_add': 107901, 'pi_1_ntt': 0.00021401, 'pi_1_mul': 0.00004611, 'pi_1_add': 0.00001426},
    # 16384: {'h_fft': 2.556, 'h_fft_ac': 0.168645, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975, 'pi_1_ntt': 0.00003884, 'pi_1_mul': 0.00001139, 'pi_1_add': 0.00000301},
    # 32768: {'h_fft': 1.267, 'h_fft_ac': 0.359552, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
}

# fft_params_m = {
#     8:   {'h_fft': 2.556, 'h_fft_ac': 0.01792, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975, 'pi_1_ntt': , 'pi_1_mul': , 'pi_1_add': },
#     16:  {'h_fft': 1.267, 'h_fft_ac': 0.032663, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
#     32:  {'h_fft': 0.623237, 'h_fft_ac': 0.065576, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
#     64:  {'h_fft': 0.301866, 'h_fft_ac': 0.163766, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
#     128: {'h_fft': 0.137299, 'h_fft_ac': 0.349483, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
#     256: {'h_fft': 0.055728, 'h_fft_ac': 0.718714, 'base_fft': 847392, 'base_mul': 150780, 'base_add': 26975},
# }


# --------------------------------------------------
# Sweep 1: range over M
# --------------------------------------------------

print("\n===== Sweep over N (M=4096, m=64) =====")

M = 4096
m = 64
N_values = [4096, 8192, 16384]

rino_h_fft = 5.134


for N in N_values:
    G, W, Gp, Wp, logG, logGp = compute_circuit_sizes(M, N, m)
    fft = fft_params_M[M]
    h_fft, h_fft_ac = fft['h_fft'], fft['h_fft_ac']
    base_fft, base_mul, base_add = fft['base_fft'], fft['base_mul'], fft['base_add']


    sub_time = runtime_sub_prover(G, W, pi_1_ntt, pi_1_mul, pi_1_add, pi_2_ntt, pi_2_mul, pi_2_add,
                                  h_fft, base_fft, base_mul, base_add)

    master_time = runtime_master_prover(Gp, Wp, pi_2_ntt, pi_2_mul, pi_2_add,
                                        h_fft_ac, base_fft, base_mul, base_add)

    verifier_time = runtime_verifier(N, M, pi_2_ntt, pi_2_mul, pi_2_add,
                                     pi_2_decrypt, pi_2_pl_mul)


    rino_time = runtime_sub_prover(G*32, W*32, pi_1_ntt, pi_1_mul, pi_1_add, pi_2_ntt, pi_2_mul, pi_2_add,
                            rino_h_fft, base_fft, base_mul, base_add)

    print(f"N={N:5d} | G={int(G):10d} | G'={int(Gp):8d} | "
          f"log2(G)={logG:5.2f} | log2(G')={logGp:5.2f} | "
          f"Sub={round(sub_time)} | Master={round(master_time)} | Verifier={round(verifier_time)} | Rino={round(rino_time)}")


