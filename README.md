# Alignment_CT
Alignment solutions in CT image processing



# Virtual Alignment Method (VAM) for projection image alignment algorithm to get clean CT image reconstructions.
1. VAM for rigid sample
2. VAM for elastic sample



# VAM Test Version Description
1. Generate sinogram using Radon Transform for Shepp-Logan phantom image
2. Apply random translation error to each column of sinogram
3. Calculate center of attenuation and use it as a fixed point
4. Apply VAM to the obtained fixed point
5. CT image reconstruction acquisition

P.S. See "fp_phantom.png" for the location of the fixed point in the reconstruction



# How to apply VAM to real samples with translation errors
1. Calculate a fixed point in the whole projection set
2. Align the projection image set to the axial level (common layer) of the specimen by adjusting the height of the fixed point
3. Apply VAM to a fixed point on the adjusted projection image set
4. Obtain ideally aligned CT reconstructions using an ideally aligned projection image set
