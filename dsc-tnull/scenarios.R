sourceDir("datamakers")
scenarios=list()

#Now, for each scenario create an element of scenarios of the following form
scenarios[[1]]=list(name="null.df4",
                    fn=null_t_datamaker,
                    args=list(
                      nsamp=1000, df=4
                      ),
                    seed=1:100)

scenarios[[2]]=list(name="null.df10",
                    fn=null_t_datamaker,
                    args=list(
                      nsamp=1000, df=10
                    ),
                    seed=1:100)

scenarios[[3]]=list(name="null.df100",
                    fn=null_t_datamaker,
                    args=list(
                      nsamp=1000, df=100
                    ),
                    seed=1:100)

scenarios[[4]]=list(name="null.df1000",
                    fn=null_t_datamaker,
                    args=list(
                      nsamp=1000, df=1000
                    ),
                    seed=1:100)