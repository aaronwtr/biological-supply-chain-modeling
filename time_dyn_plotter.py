import numpy as np
import output_handler as oh
import matplotlib.pyplot as plt


def one_d_time_series_plotter(results, sp_type):
    if sp_type == 'B':
        specie = 'bacteria'
    elif sp_type == 'L':
        specie = 'infected bacteria'
    elif sp_type == 'P':
        specie = 'phages'
    else:
        return print('ERROR: Enter specie type (B, L or P)')


    Bt = np.transpose(results)

    tot_nums = []
    for i in range(len(Bt)):
        tot_nums.append(sum(Bt[i]))

    time_series = []
    for i in range(len(tot_nums)):
        time_series.append(i)

    plt.semilogy(time_series, tot_nums)
    plt.xlabel('Time')
    plt.ylabel('Number of {}'.format(specie))
    plt.title('Time series plot of {} population change'.format(specie))
    plt.show()

    return


Bresults = oh.read_file("results/bacteria_sim_results.txt")
Lresults = oh.read_file("results/inf_bacteria_sim_results.txt")
Presults = oh.read_file("results/phages_sim_run.txt")

one_d_time_series_plotter(Bresults, 'B')
one_d_time_series_plotter(Lresults, 'L')
one_d_time_series_plotter(Presults, 'P')
