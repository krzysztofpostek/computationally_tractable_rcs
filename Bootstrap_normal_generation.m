% Bootstrap_normal_generation

N_sample_normal  = 10^5;

Error = std_deviation_common*repmat(randn(1,N_sample_normal),[N_antennas 1]) + std_deviation_intrinsic*randn(N_antennas,N_sample_normal);