feed_ratio_1 = 1/1000;
feed_ratio_2 = 3/2000;
rate_ratio = 1/0.01;
conversion_1 = 0.631;
conversion_2 = 0.794;
[Mn_1, Mw_1, PDI_1,avg_T_1, avg_DB_1] = main(rate_ratio, feed_ratio_1,conversion_1);
[Mn_2, Mw_2, PDI_2,avg_T_2, avg_DB_2] = main(rate_ratio, feed_ratio_2,conversion_2);
