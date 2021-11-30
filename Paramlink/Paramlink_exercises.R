### Exercises with Paramlink
# See documentation here: https://cran.r-project.org/web/packages/paramlink/paramlink.pdf
# Install: 
# install.packages("paramlink")

## Set up the environment
# If on the biocompace HPC server, load the R module
# 'module load R/3.6.3'
# Start an R interactive session:
# 'R'
# Load paramlink
library(paramlink)

# Set working directory
setwd("~/Projects/Seminars/ACE_Uganda/Paramlink")  # Modify this to match the location where your dowloaded files are

# Check for merlin
system("which merlin")
  # /Users/oleraj/anaconda3/bin/merlin
# If you don't see merlin but you have it installed, then run these steps:
system("echo $PATH")
# I see this:
  # /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin
# Modify this command below to prepend your path to merlin.  
# In my case merlin is installed in '/Users/oleraj/anaconda3/bin'.
# If you don't know where merlin is installed, open a terminal, 
# type 'which merlin'.  
# If you installed merlin on the biocompace HPC, you may need to initiate conda 
# for your terminal with 'eval "$(conda shell.bash hook)" '
Sys.setenv(PATH="/Users/oleraj/anaconda3/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin") 
system("which merlin")
# If you are using the biocompace HPC server, you can load merlin by loading the
# anaconda module: 'module load anaconda3-2019.10-gcc-9.3.0-7gh72na'
system("which merlin")
  # /opt/ohpc/admin/spack/0.15.0/opt/spack/linux-centos8-sandybridge/gcc-9.3.0/anaconda3-2019.10-7gh72nae2c5prph4unxinmialubphorl/bin/merlin




#### 1. Drawing simple pedigrees with Paramlink  ####
## Drawing simple pedigrees
# Autosomal recessive, homozygous
p1 <- nuclearPed(noffs=1, sex=c(1))  # sex 1 = male, 2 = female
plot(p1)
# Everyone is unaffected by default
# Make the child affected using the swapAff()
p1 <- swapAff(p1, 3)
plot(p1)  
# Add SNP markers where parents are HET and child is homozygous for one allele
m <- marker(p1, c(1, 2), c("G", "T"), 3, c("G", "G"), name="gene", chrom="1", pos="30")
p1 <- addMarker(p1, m)
plot(p1, marker=1)
# If you want to save an image of the file, run these commands below
pdf("hom_rec.pdf", height=3, width=3)
plot(p1, marker=1)
dev.off()

# Autosomal Recessive, Compound Heterozygous
p2 <- nuclearPed(noffs=2, sex=c(2))
plot(p2)
p2 <- swapAff(p2, 3)
p2
plot(p2)  
# Add markers showing each parent as HET for two separate alleles,
# proband (3) HET for both and unaffected sibling as HET for just one allele
# First allele ref G, alt T
# Second allele ref A, alt C
m1 <- marker(p2, c(1, 3), c("G", "T"), 2, c("G", "G"), 4, c("G", "G"))
m2 <- marker(p2, c(2, 3), c("A", "C"), 1, c("A", "A"), 4, c("A", "C"))
p2 <- addMarker(p2, m1)
p2 <- addMarker(p2, m2)
p2
plot(p2, marker = c(1,2))

# Save as PDF
pdf("comp_het.pdf", height=4, width=3)
plot(p2, marker= c(1,2))
dev.off()

# De novo
p3 <- nuclearPed(noffs=1, sex=c(2))
p3 <- swapAff(p3, 3)
m1 <- marker(p3, c(1, 2), c("C", "C"), 3, c("C", "T"))
p3 <- addMarker(p3, m1)
plot(p3, marker=1)

pdf("de_novo.pdf", height=3, width=3)
plot(p3, marker= 1)
dev.off()

# X_rec
p4 <- nuclearPed(noffs=1, sex=c(1))
p4 <- swapAff(p4, 3)
# Add markers.  You can either add them with the second allele blank for males,
# or just specify a single allele and it will convert it to homozygous
m1 <- marker(p4, 1, c("C", ""), 2, c("C", "T"), 3, c("T", ""))
m2 <- marker(p4, 1, c("C"), 2, c("C", "T"), 3, c("T"))
p4 <- addMarker(p4, m1)  # Adds as marker 1
p4 <- addMarker(p4, m2)  # Adds as marker 2
p4
plot(p4, marker=1)
plot(p4, marker=2)

pdf("X_rec.pdf", height=3, width=3)
plot(p4, marker= 1)
dev.off()

# See here for more complicated examples (multi-generational [addParents], multiple affected statuses, etc.): https://web.archive.org/web/20170420081914/http://folk.uio.no/magnusv/LinkageCourse/Paramlink/paramlink_intro.pdf 





#### 2. Power analysis  ####
# Adapted from https://web.archive.org/web/20170420081924/http://folk.uio.no/magnusv/LinkageCourse/Paramlink/paramlink_power.pdf

# Power estimation is done by simulating markers that are completely linked to the disease locus.
# Singlepoint LOD scores of these markers are then computed, and the maximum score is reported.

# Load a toy pedigree
x = linkdat("toy_example.ped")
plot(x)

# Set genetic model
# Basic genetic models with default complete penetrance, and dfreq set to 0.00001 :
# 1 = Autosomal dominant
# 2 = Autosomal recessive
# 3 = X-linked dominant
# 4 = X-linked recessive
# Or you can manually set parameters, e.g., 
# setModel(x, chrom="autosomal", penetrances=c(0.01, 0.9, 1.0), dfreq=0.005)
# penetrances: likelihood of getting disease with 1) 0 alleles (phenocopy), 2) 1 allele, 3) 2 alleles
x = setModel(x, model=1)  # dominant model, default parameters

set.seed(1234)

# Recall that for dominant, Max LOD ≈ 0.3(c − 1) 
# In our case there are 5 children (meioses) so Max LOD should be about 0.3 * 4 =~ 1.2
# Simulate with 500 markers
linkage.power(x, N=500) 

# Were we correct?

# Now what would our Max LOD score be if we were missing data for grandma?
  # Modify the "available" parameter to exclude #2
linkage.power(x, available=c(1, 3:8), N=500)

# As you can see, this reduces the Max LOD score to 1.0.  
# Missing information from grandma hurts the analysis.
# You can use this to determine which samples are most important to include

# What if we simulate using markers that are triallelic instead of biallelic?
linkage.power(x, available=c(1, 3:8), afreq=c(0.3, 0.3, 0.4), N=500)
# That brings up the LOD score back to the same level where it was when we had data for grandma
# This demonstrates that using multi-allelic markers can increase your power for linkage analysis 
# even without adding more samples.

# Try modifying penetrance scores
  # Penetrance scores are for zero, one, or two copies of the disease allele, respectively
x80 = setModel(x, penetrances=c(0, 0.8, 1))  # dominant with 80% penetrance for Heterozygotes (HET)
linkage.power(x80)

# Max LOD score is slightly reduced, 1.12.  
# This makes sense because we have greater uncertainty in the analysis. 
# Unaffected individuals are allowed to be carriers so it's harder to 
# trace the disease allele.

x60 = setModel(x, penetrances=c(0, 0.6, 1))  # dominant with 80% penetrance for Heterozygotes (HET)
linkage.power(x60)

# Max LOD score is reduced further, 1.058.  

# Now for one more, let's draw a recessive model from scratch 
  # Nuclear family with 8 children, 3 affected, 5 unaffected
p5 <- nuclearPed(noffs=8, sex=c(1,1,2,2,1,2,1,2))
plot(p5)
p5 <- swapAff(p5, c(5,6,8))
plot(p5)  

# Recall that the general rule for nuclear families in autosomal recessive is
# Max LOD = 0.6(a − 1) + 0.125n
# We have 3 aff and 5 unaff, so our Max LOD should be (0.6 * 2) + (0.125 * 5) = 1.2 + 0.625 = 1.825
# See if we're right:
p5 <- setModel(p5, model=2)
linkage.power(p5, N=500)

# Would adding more affected or unaffected help increase our LOD score for recessive disorder?

# What about for a dominant disorder?







#### 3. Computing LOD scores and Merlin wrapper  ####
# Exercises adapted from https://web.archive.org/web/20170408144024/http://folk.uio.no/magnusv/LinkageCourse/Paramlink/paramlink_merlin.pdf

# First check to see if paramlink recognizes merlin
x = linkdat("toy_example.ped", model=1)
merlin(x)
lod(x)
# lod() is the built-in two-point linkage analysis option in paramlink
# You should see 1.204 for both

# Now let's try a larger pedigree with multiple markers to compare
# two-point linkage (paramlink) to multipoint linkage analysis (MERLIN)
y = linkdat(dominant, model=1)
y  # Inspect the PED file data we loaded
plot(y, available=TRUE)  # available = TRUE marks in red the ones for which we have data
# What is our theoretical Max LOD for this pedigree?
linkage.power(y, N=500)  
# 3.6

# Merlin is slow with more than about 16 individuals, so let's subset to only 
# include affected individuals (and important unaffected)
y_aff = trim(y, keep="affected")
summary(y_aff)
plot(y_aff, available=TRUE)
# Theoretical Max LOD now?
linkage.power(y_aff, N=500) 
# 2.01


# Remove Mendelian errors
mendelianCheck(y_aff)   # 5 markers with Mendelian errors.  Inspect these manually in the raw data and/or remove. 
y_aff = mendelianCheck(y_aff, remove=TRUE)
mendelianCheck(y_aff)   # None, they were removed.

# Compute LOD scores with two-point linkage
single_lods = lod(y_aff)
# Now get multipoint linkage LOD scores (combine nearby markers, like a haplotype)
multi_lods = merlin(y_aff)

#dev.off()
plot(single_lods, lty=3) # lty=3 gives dashed line
par(new=T) # prepares R for a new plot in the same window
plot(multi_lods, col="blue") # the 'col' argument specifies line color

# Which method gives better results?

# Find the best peak
summary(multi_lods)
# Get all markers above a certain threshold
lod.peaks(multi_lods, threshold=2.4)
# Plot the pedigree with peak marker genotypes
plot(y, marker=313)  # Use 313 as representative

# Anything odd showing up? Any potentially problematic samples?

pdf("toy_ped_peak_marker_genotypes.pdf", width=4, height=7)
plot(y, marker=c(313, 314, 315, 316))
dev.off()

# Analyze with different penetrance values to see how it affects LOD score
y100 = setModel(y_aff, penetrances = c(0, 1, 1), dfreq=0.0001)
y80 = setModel(y_aff, penetrances = c(0, 0.8, 1), dfreq=0.0001)
y60 = setModel(y_aff, penetrances = c(0, 0.6, 1), dfreq=0.0001)

m_lods100 = merlin(y100)
m_lods80 = merlin(y80)
m_lods60 = merlin(y60)

# Plot results for all different penetrance values together
dev.off()
plot(m_lods100, col="green")
par(new=TRUE) # new plot in the same window
plot(m_lods80, col="red")
par(new=TRUE) # new plot in the same window
plot(m_lods60, col="blue")
legend("topleft", c("f1=100%", "f1=80%", "f1=60%"), lwd=2, col = c("green", "red", "blue"))
abline(h=0, col="black")




## Now a real case using PED, MAP information for markers (chr, position, etc.), and DAT file
  # These are the same files you would use with MERLIN
  # You would need a separate set of files per chromosome, run each chromosome separately
  # For this example, we'll use data for chromosome 10

# Load data
z = linkdat(ped="HSkr10.ped", dat="HSkr10.dat", map="HSkr10.map")
plot(z, available=TRUE)
# Set model: dominant with 90% penetrance for HETs
z = setModel(z, chrom="autosomal", penetrances=c(0, 0.9, 1), dfreq=0.0001)

# Remove Mendelian errors
mendelianCheck(z)   # 6 Mendelian errors
z = mendelianCheck(z, remove=TRUE)
mendelianCheck(z)   # None, they were removed.
m_lods = merlin(z)
dev.off()
plot(m_lods)
summary(m_lods)
lod.peaks(m_lods, threshold=1)

# Analyze with different penetrance values to see how it affects LOD score
z100 = setModel(z, penetrances = c(0.0001, 1, 1))
z80 = setModel(z, penetrances = c(0.0001, 0.8, 1))
z60 = setModel(z, penetrances = c(0.0001, 0.6, 1))

m_lods100 = merlin(z100)
m_lods80 = merlin(z80)
m_lods60 = merlin(z60)

# Plot results for all different penetrance values together
plot(m_lods100, col="green")
par(new=TRUE) # new plot in the same window
plot(m_lods80, col="red")
par(new=TRUE) # new plot in the same window
plot(m_lods60, col="blue")
legend("topleft", c("f1=100%", "f1=80%", "f1=60%"), lwd=2, col = c("green", "red", "blue"))
abline(h=0, col="black")


# Note that normally when you reduce the penetrance, the LOD score decreases
# For the smaller peaks that is the case.
# However for the tall peak, the LOD score increases with higher penetrance, 
# suggesting that if this is a real peak that some of the unaffected individuals 
# carry the disease allele.  Note that with 100% penetrance, the large peak is flat-lined.


# For running multiple families,
# 1) with Paramlink, you can save data for each family as an array (per chromosome)
#    and merge the data across families
# 2) OR (recommended) use Merlin natively in the terminal, which can handle multi-  
#    family PED files. http://csg.sph.umich.edu/abecasis/merlin/tour/parametric.html
