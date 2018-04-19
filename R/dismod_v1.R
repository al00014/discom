#' dismod_v1 function
#'
#' This function simulates the process behind DisMod-II (developed by WHO, available at www.who.int/evidence/dismod, for windows-PC only) and DisMod-MR (developed by IHME). Function in current forms only consider Stock and Condition (which is used in DisMod-MR), leaving out the Dead group (which is used in DisMod-II).
#' 
#' Beware that currently, this function is like a toy example, only taking in incidence data and mortality data
#' as its main input data. The reason to do so is that the author currently only came across situation where
#' incidence and mortality data are available (incidence and mortality data are already readily available in the death registry system). 
#'
#' The function in current forms only take in pops, alldeath_rate, condition_death_rate, condition_incidence_rate and duration where these arguments have an exact match in vector length with each other.
#'
#' Also, this model calculates remission with the equation: remission = 1/duration (adapted from the book by Abraham D. Flaxman, title "An Integrative Metaregression Framework for Descriptive Epidemiology", ISBN-13: 978-0295991849). 
#'
#' The real dismod software actually takes in three of the following epidemiological
#' variables: incidence, remission, case fatality (or RR), prevalence, mortality, and generates the rest via the integrating system modeling.
#'
#'
#' @param pops taking in a vector of age-gender-specific population, preferrably for a certain time (for example, in a year).
#' @param alldeath_rate taking in a vector of age-gender-specific all cause mortality rate, preferrably for a certain time (for example, in a year).
#' @param condition_death_rate taking in a vector of age-gender-specific mortality rate for a specific condition (a disease), preferrably for a certain time (for example, in a year).
#' @param condition_incidence_rate taking in a vector of age-gender-specific incidence rate for a specific condition (a disease), preferrably for a certain time (for example, in a year).
#' @param duration a vector specifying the age-gender-specific duration of the disease.
#' @param Stock an initial value for the compartmental model (for Stock population--susceptible population).
#' @param Condition an initial value for the compartmental model (for diseased population--population with condition).
#' @param start_age specifying the lower bound of x-axis.
#' @param age_range manually setting up the age label vector. In real-world data, mortality and incidence data are often specified in this way: "0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+". The default value for this parameter is to mimic the real world situation.
#' @param end_age specifying the upper bound of x-axis.
#' @param plots a boolean argument, specifying whether to visualize the compartmental framework model or not. The plots are in a ggplot fashsion.
#' @param ...
#' @return A gridExtra object, when plots=TRUE;
#' @return or a list of data containing the simulated data produced by differentiation in the modeling process, when plots=FALSE. The return data also contains prevalence estimated from the input incidence and mortality. 
#' @keywords compartmental framework, dismod
#' @export
#' @examples
#' dismod_v1(pops=population,
#'          alldeath_rate=ACD_rate,
#'          condition_death_rate=CD_rate,
#'          condition_incidence_rate=CI_rate,
#'          duration=CDuration,
#'          Stock=10000,
#'          Condition=5,
#'          start_age=0,
#'          age_range=c(0,1,seq(5,85,by=5)),
#'          end_age=85)         
#' 
#' 
#'
#'
							   
dismod_v1<-function(pops,
                    alldeath_rate,
                    condition_death_rate,
                    condition_incidence_rate,
                    duration,
                    Stock=10000,
                    Condition=1,
                    start_age=0,
                    age_range=c(0,1,seq(5,85,by=5)),
                    end_age=85,
                    plots=TRUE){
        message('Assuming duration=1/hr !!')
        if(length(age_range)!=length(pops)){
                stop('Please ensure age_range complies with input data!')
        }
        if(length(age_range)!=length(alldeath_rate)){
                stop('Please ensure age_range complies with input data!')
        }
        if(length(age_range)!=length(condition_death_rate)){
                stop('Please ensure age_range complies with input data!')
        }
        if(length(age_range)!=length(condition_incidence_rate)){
                stop('Please ensure age_range complies with input data!')
        }
        
        inits<-c(S=Stock,C=Condition)
        time_var=seq(start_age,end_age,by=(end_age+length(start_age))/(2*length(start_age:end_age)))
        #require(deSolve)
        equation<-function(time_var,parameters,state){
                with(as.list(c(state,parameters)),{
                        dS=-(hi+hm)*S+hr*C
                        dC=hi*S-(hr+hm+hf)*C
                        return(list(c(dS,dC)))
                })
        }
        #parameters=c(hi=0.01,
        #             hm=0.4,
        #             hr=0.3,
        #             hf=0.05)
        
        res <- vector(length(pops),mode="list")
        for(i in seq_along(pops)){
                hmall<-alldeath_rate[i]
                hm<-condition_death_rate[i]
                hi<-condition_incidence_rate[i]
                hf<-(hmall-hm)*(Stock+Condition)/Condition
                hr<-1/duration[i]
                res[[i]]<-deSolve::ode(y=inits,times=time_var,func = equation,
                              parms = c(hi=hi,
                                        hm=hm,
                                        hr=hr,
                                        hf=hf)
                              
                              )
        }
        #cout_index<-1
        new_df<-data.frame()
        for(i in 1:length(res)){
                new_df<-rbind(new_df,
                              res[[i]][res[[i]][,'time']==age_range[i],])
                #cout_index<-cout_index+1
        }
        
        colnames(new_df)<-c('age','S','C')
        new_df$prevalence_adj<-new_df$C/(new_df$S+new_df$C)
        
        title1 <- bquote("Incidence rate")
        #subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
        #))
        plot1_data<-data.frame(age=age_range,
                               y=condition_incidence_rate)
        
        plot1<-ggplot2::ggplot(plot1_data)+
                ggplot2::ggtitle(bquote(atop(bold(.(title1))#,atop(bold(.(subtit1)))
                )
                )
                )+
                ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                ggplot2::ylab(label="Rates")+
                ggplot2::xlab(label="Age (years)")#+
                #ggplot2::ylim(c(0,y_lim))
        title2 <- bquote("Remission rate")
        #subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
        #))
        plot2_data<-data.frame(age=age_range,
                               y=1/duration)
        plot2<-ggplot2::ggplot(plot2_data)+
                ggplot2::ggtitle(bquote(atop(bold(.(title2))#,atop(bold(.(subtit1)))
                )
                )
                )+
                ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                ggplot2::ylab(label="Rates")+
                ggplot2::xlab(label="Age (years)")
        ##
        title3 <- bquote("Mortality rate")
        #subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
        #))
        plot3_data<-data.frame(age=age_range,
                               y=condition_death_rate)
        plot3<-ggplot2::ggplot(plot3_data)+
                ggplot2::ggtitle(bquote(atop(bold(.(title3))#,atop(bold(.(subtit1)))
                )
                )
                )+
                ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                ggplot2::ylab(label="Rates")+
                ggplot2::xlab(label="Age (years)")
        
        ##
        title4 <- bquote("Fatality rate")
        #subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
        #))
        plot4_data<-data.frame(age=age_range,
                               y=(alldeath_rate-condition_death_rate)*(Stock+Condition)/Condition)
        plot4<-ggplot2::ggplot(plot4_data)+
                ggplot2::ggtitle(bquote(atop(bold(.(title4))#,atop(bold(.(subtit1)))
                )
                )
                )+
                ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                ggplot2::ylab(label="Rates")+
                ggplot2::xlab(label="Age (years)")
        
        
        title5 <- bquote("Two compartmental model")
        plot5<-ggplot2::ggplot(new_df)+
                ggplot2::ggtitle(bquote(atop(bold(.(title5))#,atop(bold(.(subtit3)))
                )
                )
                )+
                ggplot2::geom_line(ggplot2::aes(y=prevalence_adj,x=age))+#+theme_bw()
                ggplot2::ylab(label="Prevalence rate")+
                ggplot2::xlab(label="Age (years)")    
        
        #require(gridExtra)
        
        if(plots==TRUE){
                final_plots<-gridExtra::grid.arrange(plot1,plot2,plot5,
                                          plot4,plot3,
                                          layout_matrix=matrix(c(1,2,3,4,5,3),byrow = TRUE,ncol=3))
                print(final_plots)
                return(list(adjusted_data=new_df))
        } else{
                message('No plots are shown!')
                return(list(adjusted_data=new_df))
        }
        
        #return(list(plot_ls=final_plots,
        #            adjusted_data=new_df))
}