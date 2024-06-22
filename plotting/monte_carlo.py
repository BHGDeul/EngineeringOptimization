import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

def create_boxplot(data,objvals):

    fig, ax = plt.subplots(3, 2, figsize=(10,12))

    sns.boxplot(data=data[0], ax=ax[0,0])
    ax[0,0].set_title('Projected area (m^2')

    sns.boxplot(data=data[1], ax=ax[0,1])
    ax[0,1].set_title('Gamma_in')

    sns.boxplot(data=data[2], ax=ax[1,0])
    ax[1,0].set_title('Gamma_out')

    sns.boxplot(data=data[3], ax=ax[1,1])
    ax[1,1].set_title('Elevation in (rad)')

    sns.boxplot(data=data[4], ax=ax[2,0])
    ax[2,0].set_title('Elevation out (rad)')

    sns.boxplot(data=objvals, ax=ax[2,1])
    ax[2,1].set_title('Objective function value')

    fig.suptitle('Monte-Carlo analysis of minimizers')

    plt.show()
    plt.savefig('output/plots/monte_carlo_cobyla.png')