# Loading anaconda
module load anaconda3/2022.10

env_dir="codes/GCF_pipeline"

# Checking if conda environment has already been created, and creating if it hasn't

if [ -d "${env_dir}" ]; then
	echo "Conda environment has been previously created."
else
	echo "Conda environment has not been previously created. Creating."
	conda env create -p codes/GCF_pipeline -f codes/GCF_pipeline.yml
fi
