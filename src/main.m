%% MRI Fuzzy Segmentation of Brain using Neighborhood Extraction
%
% Paper Authors: Shan Shen, William Sandham, et al.
%
% Implementation by: Haseeb, Bhavishey, Dishant
%
%
%% 1) Formulation
% 
% Precise segmentation of the brain is important for detection of tumors,
% edema, and necrotic tissues so that proper diagnosis and treatment can be 
% provided. However, segmentation techniques such as thresholding, which rely
% heavily on contrast resolution, are prone to limitations due to the 
% complex distribution of tissue intensities in MRI images. 
% The large number of intensities and the presence of noise makes it difficult 
% to determinine appropriate thresholds.
%
% As a result of these limitations, thresholding must be combined with other 
% segmentation methods such as fuzzy-c-means (FCM) clustering or 
% expectation-maximization (EM). However, both FCM and EM are themselves prone
% to the effects of noise. Traditional denoising methods which perform Statistical 
% estimation of noise, as being Gaussian or Laplacian for example, are
% invalid during MRI segmentation. In fact, the process of segmentation itself
% introduces further noise.
%
% A second problem associated with FCM is the lack of consideration for spatial 
% dependence. Pixels are associated with a cluster based solely on pixel 
% intensity values. Distance of a pixel to a cluster is not featured in the membership, 
% and therefore, noise has a paramount effect on the result of segmentation.
% Hence, clustering that provides robust and 
% consistent segmentation in the presence of unknown noise is desired.
% 
% 
%
%
%% 2) Proposed Solution
%
% The algorithm proposed in the paper, termed "improved fuzzy c-means" (IFCM), reduces 
% the effect of noise by considering a neighborhood of pixel intensities
% and distances to a centroid. This is in contrast to FCM which considers only 
% individual pixel intensities to a centroid. This algorithm improves two key
% factors in Neighborhood Attraction:
%
% * feature attraction: differences between neighboring pixel intensities
% * distance attraction: relative locations of neighboring pixel
% 
% The dataset $\mathbf{X}$ is a _pxN_ matrix where _p_ is the dimension of each
% $\mathbf{x}_{j}$ vector and _N_ is the size of the image. For the basic FCM case,
% _p_ is taken to be 1 such that each $\mathbf{x}_{j}$ vector is _1x1_ and holds only the 
% intensity of a pixel. $\mathbf{X}$ is thus constructed by column-stacking the input
% image.
%
% The membership function of vector $\mathbf{x}_j$ to ${i}^{th}$ cluster is given by:    
% 
% $$\mathbf{u}_{ij} =\frac{1}{ \sum_{k=1}^{C}(\frac{\mathrm{d (x_j{}, v_i{})}
% }{\mathrm{d} (x_j{}, v_k{}) })^{2/(m-1)} }$$
%
% The ${i}^{th}$ cluster feature center is:
%
% $$\mathbf{v}_i = \frac{ \sum_{j=1}^{N}(u_{ij})^{m} \mathbf{x}_j }{ 
% \sum_{j=1}^{N}(u_{ij})^{m} }$$
%
% The primary differentiator between IFCM and FCM is IFCM's consideration of 
% Neighborhood Attraction in the degree of fuzziness as follows:
%
% $$\mathrm{d^{2}} (\mathbf{x}_j , \mathbf{v}_j) = ||\mathbf{x}_j - \mathbf{v}_j
% || (1 - \lambda \mathbf{H}_{ij} - \xi \mathbf{F}_{ij})$$
%
% The matrix $\mathbf{H}_{ij}$ is constructed using the intensity difference between
% pixel j and its neighbour k:
%
% $$\mathbf{H}_{ij} = \frac{ \sum_{k=1}^{S} \mathbf{u}_{ik} || \mathbf{x}_j - 
% \mathbf{x}_k ||}{ \sum_{k=1}^{S} || \mathbf{x}_j - \mathbf{x}_k || }$$
%
% Likewise, the matrix $\mathbf{F}_{ij}$ is constructed using the spatial difference,
% $q_{jk}$, which is governed by the relative location of pixel k to pixel j within the
% neighborhood.
%
% $$\mathbf{F}_{ij} = \frac{ \sum_{k=1}^{S} \mathbf{u}_{ik}^2 \mathbf{q}_{jk}^2 }
% { \sum_{k=1}^{S} \mathbf{q}_{jk}^2 }$$
%
% Finally, clusters are decided iteratively updating U and V to minimize the following cost function:
%
% $$\mathit{argmin} \sum_{j=1}^{N} \sum_{i=1}^C \mathbf{u}_{ij}^2 d^2(\mathbf{x}_j,\mathbf{v}_i)$$
%
% The stop condition is decided by a degree of convergence such that 
% $| U^{(l+1)}-U^{(l)} | \leq \varepsilon$ where $l$ represents the loop iteration and 
% $\varepsilon$ is the convergence parameter.
%
% Based on the research paper, the results presented here use $\lambda = 0.6038 $ and 
% $\varepsilon = 0.6097$
%
%% 3) Data Sources
%
% A brain MRI with a tumour and noticable salt-pepper noise is acquired from 
% Figure3 of <http://www.ajnr.org/content/27/3/475/tab-figures-data>. This image
% is purely grayscale.
%
% This image is then subjected to additive noise with a tunable parameter, alpha
% to produce 2 noisy copies.
clc;
clear;
      
f1 = im2double(imread('brain-tumour-mri.gif'));
f1=f1(135:250,215:330);    
subplot(1,3,1), imshow(f1);
title('f1 (original)');

alpha = 0.4;
noise = rand(size(f1));
f2 = f1 + alpha*noise;
subplot(1,3,2), imshow(f2);
title('f2 (noisy)');

alpha = 0.7;
f3 = f1 + alpha*noise;
subplot(1,3,3), imshow(f3);
title('f3 (very noisy)');
    
%% 4) Solution
% 
% The membership _U_ and the centroids _V_, are computed iteratively
% through FCM and IFCM on each of the test images: f1, f2, and f3 respectively.
%
% Details and source code of the implementation of FCM and IFCM are available in the 
% 'custom source files' section. The IFCM follows the same order of steps as presented in 
% the paper. In summary, for IFCM, a full round of FCM is executed to determine the
% first membership matrix U for IFCM. Then the degree of fuzziness is computed with
% the matrices $H_{ij}$ and $F_{ij}$ along with parameters $\lambda$ and $\varepsilon$.
% These matrices are computed using neighborhood kernel operations.
%
%


% Segment dataset for image f1
X1 = f1(:);  % Load Dataset
[U1_fcm, V1_fcm] = fcm(X1, 3, 2);     % Execute FCM on f1
maxU = max(U1_fcm);
index1_fcm = find(U1_fcm(1,:) == maxU);
index2_fcm = find(U1_fcm(2, :) == maxU);
index3_fcm = find(U1_fcm(3, :) == maxU);

[U1_ifcm, V1_ifcm] = ifcm(X1, 3, 2);  % Execute IFCM on f1
maxU = max(U1_ifcm);
index1_ifcm = find(U1_ifcm(1,:) == maxU);
index2_ifcm = find(U1_ifcm(2, :) == maxU);
index3_ifcm = find(U1_ifcm(3, :) == maxU);

% Segment dataset for image f2
X2 = f2(:);  % Load Dataset
[U2_fcm, V2_fcm] = fcm(X2, 3, 2);     % Execute FCM on f2
maxU = max(U2_fcm);
index1_fcm2 = find(U2_fcm(1,:) == maxU);
index2_fcm2 = find(U2_fcm(2, :) == maxU);
index3_fcm2 = find(U2_fcm(3, :) == maxU);

[U2_ifcm, V2_ifcm] = ifcm(X2, 3, 2);  % Execute IFCM on f2
maxU = max(U2_ifcm);
index1_ifcm2 = find(U2_ifcm(1,:) == maxU);
index2_ifcm2 = find(U2_ifcm(2, :) == maxU);
index3_ifcm2 = find(U2_ifcm(3, :) == maxU);

% Segment dataset for image f3
X3 = f3(:);   % Load Dataset
[U3_fcm, V3_fcm] = fcm(X3, 3, 2);     % Execute FCM on f3
maxU = max(U3_fcm);
index1_fcm3 = find(U3_fcm(1,:) == maxU);
index2_fcm3 = find(U3_fcm(2, :) == maxU);
index3_fcm3 = find(U3_fcm(3, :) == maxU);

[U3_ifcm, V3_ifcm] = ifcm(X3, 3, 2);  % Execute FCM on f3
maxU = max(U3_ifcm);
index1_ifcm3 = find(U3_ifcm(1,:) == maxU);
index2_ifcm3 = find(U3_ifcm(2, :) == maxU);
index3_ifcm3 = find(U3_ifcm(3, :) == maxU);

%% 5) Visualization of Results
%

f1_fcm(1:length(X1))=0;
f1_fcm(index1_fcm)=1;
f1_fcm(index2_fcm)=0.5;
f1_fcm(index3_fcm)=0.2;
f1_fcm_res = reshape(f1_fcm, 116, 116);

f1_ifcm(1:length(X1))=0;
f1_ifcm(index1_ifcm)=1;
f1_ifcm(index2_ifcm)=0.5;
f1_ifcm(index3_ifcm)=0.2;
f1_ifcm_res = reshape(f1_ifcm, 116, 116);

figure;                             % New figure window for results
subplot(3,2,1), imshow(f1_fcm_res);
title('f1 FCM');

subplot(3,2,2), imshow(f1_ifcm_res);
title('f1 IFCM');

f2_fcm(1:length(X2))=0;
f2_fcm(index1_fcm2)=1;
f2_fcm(index2_fcm2)=0.5;
f2_fcm(index3_fcm2)=0.2;
f2_fcm_res = reshape(f2_fcm, 116, 116);

f2_ifcm(1:length(X2))=0;
f2_ifcm(index1_ifcm2)=1;
f2_ifcm(index2_ifcm2)=0.5;
f2_ifcm(index3_ifcm2)=0.2;
f2_ifcm_res = reshape(f2_ifcm, 116, 116);


subplot(3,2,3), imshow(f2_fcm_res);
title('f2 FCM');

subplot(3,2,4), imshow(f2_ifcm_res);
title('f2 IFCM');

f3_fcm(1:length(X3))=0;
f3_fcm(index1_fcm3)=1;
f3_fcm(index2_fcm3)=0.5;
f3_fcm(index3_fcm3)=0.2;
f3_fcm_res = reshape(f3_fcm, 116, 116);

f3_ifcm(1:length(X3))=0;
f3_ifcm(index1_ifcm3)=1;
f3_ifcm(index2_ifcm3)=0.5;
f3_ifcm(index3_ifcm3)=0.2;
f3_ifcm_res = reshape(f3_ifcm, 116, 116);

subplot(3,2,5), imshow(f3_fcm_res);
title('f3 FCM');

subplot(3,2,6), imshow(f3_ifcm_res);
title('f3 IFCM');


%% 6) Analysis and Conclusion
% 
%
%% 1) Have you been able to reproduce results reported in paper?
%  
% Yes, to some extent, with a few compromises. 
% The program requires high RAM to execute and crashed with images of size 256x256.
% Our tests were run on an MRI grayscale image of size 116x116. 
% The Neighborhood Attraction implementation also uses kernels of size 3x3 as opposed to
% 9x9 as in the paper. The expansion of these kernels to a larger size is trivial but 
% we have concerns about the program completing.
%
% The major claims of the paper have been verified through inspection. Unfortunately,
% overlap metric (comparison score) used in the paper does not clearly define what was 
% used as the 'ground truth', and has hence not been reproduced here. 
% Overall, however, the IFCM segmentation is far more 
% resistent to noise than the FCM alone.
%
% 
%% 2) Did algorithm behave in a predictable way (as described by authors)?
%  
% The IFCM algorithm is definitely more robust than FCM, as suggested by the
% authors. However, a segmented result free of noise is not guaranteed by either
% FCM or IFCM. Our results show that for a particular choice of initial membership
% the result of IFCM is clearly superior to that of FCM in all noise cases.
% Trials with FCM demonstrate that FCM usually introduces more noise but in some runs,
% IFCM and FCM have similar results, particularly in the 'very noisy' case.
%
% In general, IFCM results in a final segmentation with almost 
% no salt and pepper noise for the original and noisy image. The very noisy image
% however still has a great degree of noise in the final result for both FCM and
% IFCM. 
%
% A major observation is that the quality of segmentation is heavily influenced by the 
% initial membership matrix. This is not as evident with IFCM, possibly due to
% the fact that it goes through further iterations to determine the final membership,
% including a full round of FCM. However, the choice of membership can cause the final
% segmented image to be either extremely clear or completely unresolvable. Currently,
% the initial membership is generated randomly as the paper did not clearly state
% what was used for their tests.
%
%
%% 3) Do your own conclusions support those made by the authors?
%
% Our conlusion does support those of the authors on the condition that a suitable
% first approximation of the membership is made. A major observation is that
% the IFCM performs better segmentation than FCM for a wider range of first approximations of U.
% In general, the result of IFCM is also free of salt and pepper noise for the first two
% test cases. However, the result of FCM and IFCM tend to be very similar for the 
% 'very noisy' test case when the first membership is constructed using 70% Gaussian 
% noise.
% 
% Secondly, the choice of $\lambda$ and $\varepsilon$ do indeed make a huge difference
% in the quality of segmentation using IFCM. A poor choice of these constants can 
% cause obscure and noisy segmentation. We weren't able to pin-point the best choice
% of these values.
%  
%
%
%% 4) What are the drawbacks (if any) of the proposed solution? 
% 
% A major drawback of IFCM, and FCM for that matter, is that it requires an
% an appropriate initial approximation of the membership function. 
% Neighborhood operations are also very expensive and occupy the bulk of the 
% execution time.
% Finally, the choice of $\lambda$ and $\varepsilon$ have a huge impact on the result
% and require advanced methods such as ANN to find optimal values for each image.
%
%% 7) Custom source files
% 
% * FCM <fcm.html>
% * IFCM <ifcm.html>
% * DIFF INTENSITY <kernel_diffIntensity.html>
% * EUCLIDEAN DISTANCE <kernel_distEuclidean.html>
%
% 
%
%