import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from GMM_DP import GaussianMixtureDP
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)


def Kmeans_clu_method(cov_matrix, sample_list, pca_components, n_clusters, output_dir, save_figures, seed, gmm_dp):
    cov_var = np.var(cov_matrix, axis=0)
    pca = PCA(n_components=pca_components)
    X_new = pca.fit_transform(cov_matrix.T)
    if save_figures:
        plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=cov_var, cmap="afmhot_r")
        plt.colorbar()
        plt.savefig(os.path.join(output_dir, "cell_variance.png"))
        plt.close()
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed).fit(X_new)
    if save_figures:
        scatter = plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=kmeans.labels_)
        plt.legend(handles=scatter.legend_elements()[0], labels=[i for i in range(n_clusters)])
        plt.savefig(os.path.join(output_dir, "cell_clustering.png"))
        plt.close()
    
    cluster_variance_list = []
    for i in range(n_clusters):
        cov_matrix_c = cov_matrix[:, kmeans.labels_==i]
        cov_profile = np.median(cov_matrix_c, axis=1)
        cluster_variance_list.append(np.var(cov_profile))
        logger.info("Cluster %s: samples %s variance %s;"%(i, cov_matrix_c.shape[1], cluster_variance_list[-1]))
    idx = cluster_variance_list.index(min(cluster_variance_list))
    logger.info("Identify cluster%s as the normal cells."%idx)
    
    normal_profile = np.median(cov_matrix[:, kmeans.labels_==idx], axis=1)
    cov_matrix -= np.tile(normal_profile, (cov_matrix.shape[1],1)).T
    normal_sample_list = list(np.array(sample_list)[kmeans.labels_==idx])

    return cov_matrix, normal_sample_list
   

def Ward_clu_method(cov_matrix, sample_list, pca_components, n_clusters, output_dir, save_figures, seed, gmm_dp):
    cov_var = np.var(cov_matrix, axis=0)
    pca = PCA(n_components=pca_components)
    X_new = pca.fit_transform(cov_matrix.T)
    if save_figures:
        plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=cov_var, cmap="afmhot_r")
        plt.colorbar()
        plt.savefig(os.path.join(output_dir, "cell_variance.png"))
        plt.close()
    
    ward = AgglomerativeClustering(linkage="ward", n_clusters=n_clusters).fit(X_new)
    if save_figures:
        scatter = plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=ward.labels_)
        plt.legend(handles=scatter.legend_elements()[0], labels=[i for i in range(n_clusters)])
        plt.savefig(os.path.join(output_dir, "cell_clustering.png"))
        plt.close()
    
    cluster_variance_list = []
    for i in range(n_clusters):
        cov_matrix_c = cov_matrix[:, ward.labels_==i]
        cov_profile = np.median(cov_matrix_c, axis=1)
        if gmm_dp:
            gm = GaussianMixtureDP(n_components=3).fit(cov_profile)
            weight = list(gm.weights_)
            mean = list(abs(gm.means_))
            cluster_variance_list.append(gm.variances_[mean.index(min(mean))])
        else:
            gm = GaussianMixture(n_components=3, random_state=seed).fit(cov_profile.reshape(-1, 1))
            weight = list(gm.weights_)
            mean = list(abs(gm.means_.reshape(-1)))
            cluster_variance_list.append(gm.covariances_[mean.index(min(mean)), 0, 0])
        logger.info("Cluster %s: samples %s variance %s;"%(i, cov_matrix_c.shape[1], cluster_variance_list[-1]))
    idx = cluster_variance_list.index(min(cluster_variance_list))
    
    logger.info("Identify cluster%s as the normal cells."%idx)
    normal_profile = np.median(cov_matrix[:, ward.labels_==idx], axis=1)
    cov_matrix -= np.tile(normal_profile, (cov_matrix.shape[1],1)).T
    normal_sample_list = list(np.array(sample_list)[ward.labels_==idx])


    logger.info("Identify %s normal cells."%len(normal_sample_list))
    
    if len(normal_sample_list) < 5:
        raise ValueError("The identified normal cells is less than five cells.")

    return cov_matrix, normal_sample_list



def GMM_method(cov_matrix, sample_list, pca_components, n_clusters, output_dir, save_figures, seed, gmm_dp):
    cov_var = np.var(cov_matrix, axis=0)
    pca = PCA(n_components=pca_components)
    X_new = pca.fit_transform(cov_matrix.T)
    if save_figures:
        plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=cov_var, cmap="afmhot_r")
        plt.colorbar()
        plt.savefig(os.path.join(output_dir, "cell_variance.png"))
        plt.close()
   
    GMM_result = []
    for i in range(cov_matrix.shape[1]):
        gene_profile = cov_matrix[:, i]
        if gmm_dp:
            gm = GaussianMixtureDP(n_components=3).fit(gene_profile)
            ratio = np.sum(gm.weights_[abs(gm.means_)<0.05])
        else:
            gm = GaussianMixture(n_components=3, random_state=seed).fit(gene_profile.reshape(-1, 1))
            ratio = np.sum(gm.weights_[abs(gm.means_.reshape((-1)))<0.05])
        GMM_result.append(ratio)
    GMM_result = np.array(GMM_result)
   
    normal_profile = np.median(cov_matrix[:, GMM_result>0.8], axis=1)
    cov_matrix -= np.tile(normal_profile, (cov_matrix.shape[1], 1)).T
    normal_sample_list = list(np.array(sample_list)[GMM_result>0.8])

    if save_figures:
        scatter = plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=GMM_result, cmap="afmhot_r")
        plt.colorbar()
        plt.savefig(os.path.join(output_dir, "cell_gmm.png"))
        plt.close()
   
    logger.info("Identify %s normal cells."%len(normal_sample_list))
    
    if len(normal_sample_list) < 5:
        raise ValueError("The identified normal cells is less than five cells.")
    
    return cov_matrix, normal_sample_list



def synthetic_method(cov_matrix, sample_list, pca_components, n_clusters, output_dir, save_figures, seed, gmm_dp):
    
    logger.info("Synthesis normal cells.")

    cov_var = np.var(cov_matrix, axis=0)
    pca = PCA(n_components=pca_components)
    X_new = pca.fit_transform(cov_matrix.T)
    if save_figures:
        plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=cov_var, cmap="afmhot_r")
        plt.colorbar()
        plt.savefig(os.path.join(output_dir, "cell_variance.png"))
        plt.close()
    
    ward = AgglomerativeClustering(linkage="ward", n_clusters=n_clusters).fit(X_new)
    if save_figures:
        scatter = plt.scatter(X_new[:, 0], X_new[:, 1], marker='o', c=ward.labels_)
        plt.legend(handles=scatter.legend_elements()[0], labels=[i for i in range(n_clusters)])
        plt.savefig(os.path.join(output_dir, "cell_clustering.png"))
        plt.close()
    
    cluster_variance_list = []
    for i in range(n_clusters):
        cov_matrix_c = cov_matrix[:, ward.labels_==i]
        var_list = np.std(cov_matrix_c, axis=1)
        normal_profile = []
        for var in var_list:
            normal_profile.append([np.random.normal(loc=0, scale=var, size=1)[0]])
        cov_matrix[:, ward.labels_==i] -= np.tile(np.array(normal_profile), cov_matrix_c.shape[1])
        
    normal_sample_list = []

    return cov_matrix, normal_sample_list


