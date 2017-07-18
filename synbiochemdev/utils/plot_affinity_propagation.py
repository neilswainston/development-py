'''
=================================================
Demo of affinity propagation clustering algorithm
=================================================

Reference:
Brendan J. Frey and Delbert Dueck, 'Clustering by Passing Messages
Between Data Points', Science Feb. 2007

'''
from itertools import cycle
import random

from sklearn.cluster import AffinityPropagation
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, \
    completeness_score, homogeneity_score, silhouette_score, v_measure_score
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def get_data():
    '''Gets data.'''
    centers = [[random.random() * 2 - 1,
                random.random() * 2 - 1,
                random.random() * 2 - 1]
               for _ in range(random.randint(2, 5))]

    print len(centers)
    return make_blobs(n_samples=len(centers) * 100, centers=centers,
                      cluster_std=0.1, random_state=0)


def cluster(x_data, labels_true):
    '''Cluster.'''
    aff_prop = AffinityPropagation().fit(x_data)
    cluster_centers_indices = aff_prop.cluster_centers_indices_
    labels = aff_prop.labels_

    n_clusters = len(cluster_centers_indices)

    print 'Estimated number of clusters: %d' % n_clusters
    print 'Homogeneity: %0.3f' % homogeneity_score(labels_true, labels)
    print 'Completeness: %0.3f' % completeness_score(labels_true, labels)
    print 'V-measure: %0.3f' % v_measure_score(labels_true, labels)
    print 'Adjust Rand Index: %0.3f' % adjusted_rand_score(labels_true,
                                                           labels)
    print 'Adjust Mutual Info: %0.3f' % adjusted_mutual_info_score(labels_true,
                                                                   labels)
    print 'Silhouette Coeff: %0.3f' % silhouette_score(x_data, labels,
                                                       metric='sqeuclidean')

    return n_clusters, labels, cluster_centers_indices


def plot(x_data, n_clusters, labels, cluster_centers_indices):
    '''Plot.'''
    plt.close('all')
    fig = plt.figure(1)
    plt.clf()
    axes = fig.gca(projection='3d')

    colors = cycle('bgrcmyk')

    for k, col in zip(range(n_clusters), colors):
        class_members = labels == k
        cluster_center = x_data[cluster_centers_indices[k]]
        axes.plot(x_data[class_members, 0],
                  x_data[class_members, 1],
                  x_data[class_members, 2],
                  col + '.')

        for x_val in x_data[class_members]:
            plt.plot([cluster_center[0], x_val[0]],
                     [cluster_center[1], x_val[1]],
                     [cluster_center[2], x_val[2]],
                     col)

        axes.plot([cluster_center[0]],
                  [cluster_center[1]],
                  [cluster_center[2]],
                  'o',
                  markerfacecolor=col,
                  markeredgecolor='k', markersize=10)

    plt.title('Estimated number of clusters: %d' % n_clusters)
    plt.show()


def main():
    '''main method'''
    x_data, labels_true = get_data()
    n_clusters, labels, cluster_centers_indices = cluster(x_data, labels_true)
    plot(x_data, n_clusters, labels, cluster_centers_indices)


if __name__ == '__main__':
    main()
