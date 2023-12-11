#!/usr/bin/env bash
guppy_basecaller --input_path input_path --save_path save_path -c dna_r9.4.1_450bps_hac.cfg --device "cuda:0" --chunk_size 3000 --gpu_runners_per_device 8 --chunks_per_runner 512 --barcode_kits EXP-NBD104
