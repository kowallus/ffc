#Listeria monocytogenes bacterial genomes were downloaded from US National Center for Biotechnology Information:
#https://www.ncbi.nlm.nih.gov/pathogens

echo "1. compress genome"
ffc -f -i GCA_000585755.1_Lm1823_genomic.fna -o GCA_000585755.1_Lm1823_genomic.ffc

echo
echo "2. decompress to 'out' folder"
mkdir out
ffc -d -i GCA_000585755.1_Lm1823_genomic.ffc -o out/GCA_000585755.1_Lm1823_genomic.fna

echo
echo "3. validation"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic.fna

echo
echo "4. decompress to 'out' folder using python decompressor"
python3 ffc_decompressor.py GCA_000585755.1_Lm1823_genomic.ffc out/GCA_000585755.1_Lm1823_genomic_py.fna

echo
echo "5. validation"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic_py.fna

echo
echo "6. cleanup"
rm -R out

