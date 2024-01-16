module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/MeioticDrive2022/nextflow/main.nf \
	--fasta_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome/OMAstalkie/ \
	--gtf_dir /mnt/parscratch/users/bop20pp/Drive/ref_genome/OMAstalkie/ \
	--cellranger /users/bop20pp/software/cellranger-7.2.0/bin/cellranger \
	--metadata /users/bop20pp/personal_git/MeioticDrive2022/nextflow/metadata_mt.csv \
	--read_dir /mnt/parscratch/users/bop20pp/Drive/reads \
	--cellcycle_markers /users/bop20pp/personal_git/MeioticDrive2022/nextflow/data/cell_cycle_markers.csv \
	--celltype_markers /users/bop20pp/personal_git/MeioticDrive2022/nextflow/data/markers_elife_plus_extras.csv \
	-resume \
	-with-dag flowchat.png \
	-with-trace




#--quickgtf /mnt/parscratch/users/bop20pp/Drive/ref_genome/stalkie_mt/stalkie_mt.gtf \
        #--alevin_genome /users/bop20pp/data/Genomes/Stalkie_Dec_2023/alevin_input \
