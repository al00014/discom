#' compartmental_framework_single function
#'
#' This function estimates a single compartmental system with only Stock, birth_rate (flowing in rate) and 
#' mortality_rate (flowing out rate). This function is sort of a toy example for setting up the background of the package.
#'
#' @param birth_rate a parameter for model, can be a number or a range of numbers, must NOT be a negative value.
#' @param mortality_rate a parameter for model, can be a number or a range of numbers, must NOT be a negative value.
#' @param constant a bool parameter indicating whether the input parameters are a range of numbers or just one single number.
#' @param Stock an initial value for the compartmental model.
#' @param age specifying the extend of x-axis.
#' @param ...
#' @return A gridExtra object, and also print out the compartmental model in ggplot style.
#' @keywords compartmental framework
#' @export
#' @examples
#' Example 1: (taking only one number for parameters)
#' compartmental_framework_single(birth_rate=0.01,
#'                                mortality_rate=0.015,
#'                                Stock=1,
#'                                age=85)
#' 
#' 
#' 
#' Example 2: (taking a range of numbers for parameters) 
#'compartmental_framework_single(birth_rate=seq(0.05,0.01,length.out = length(seq(1,85,by=0.5))),
#'                               mortality_rate=seq(0.01,0.015,length.out = length(seq(1,85,by=0.5))),
#'                               constant=FALSE,
#'                               Stock=1,
#'                               age=85)
#'
#'
							   
compartmental_framework_single<-function(birth_rate,
                                         mortality_rate,
                                         constant=TRUE,
                                         Stock=1,
                                         age=85){
        
        ##### setting up plot themes
        require(ggplot2)
        mytheme4 <- theme_bw() +
                theme(text=element_text(colour="black")) +
                theme(panel.grid = element_line(colour = "white")) +
                theme(panel.background = element_rect(fill = "white"))
        theme_set(mytheme4)
        
        #### setting up necessary arguments
        time_var=seq(1,age,by=age/(2*length(1:age)))
        inits<-c(S=Stock)
        equation<-function(time_var,parameters,state){
                with(as.list(c(state,parameters)),{
                        dS=hb*S-hm*S
                        return(list(c(dS)))
                })
        }
        
        require(deSolve)
        if(constant==TRUE){
                
                parameters=c(hb=birth_rate,hm=mortality_rate)
                
                
                #equation<-function(time_var,parameters,state){
                #        with(as.list(c(state,parameters)),{
                #                dS=hb*S-hm*S
                #                return(list(c(dS)))
                #        })
                #}
                outs<-deSolve::ode(y=inits,times=time_var,func = equation,parms = parameters)
                
                out.df<-as.data.frame(outs)
                
                plot1_data<-data.frame(age=time_var,y=rep(birth_rate,times=length(time_var)))
                title1 <- bquote("Birth rate")
                subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
                ))
                
                plot2_data<-data.frame(age=time_var,y=rep(mortality_rate,times=length(time_var)))
                title2 <- bquote("Mortality rate")
                subtit2 <- bquote(list(hm==.(parameters[2])#,~gamma==.(parameters[2])
                ))
                
                
                title3 <- bquote("Single compartmental model")
                subtit3 <- bquote(list(hb==.(parameters[1]),~hm==.(parameters[2])
                ))
                if(birth_rate>mortality_rate){
                        y_lim=birth_rate
                } else {
                        y_lim=mortality_rate
                }
                plot1<-ggplot2::ggplot(plot1_data)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title1)),atop(bold(.(subtit1))))))+
                        ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                        ggplot2::ylab(label="Rates")+
                        ggplot2::xlab(label="Age (years)")+
                        ggplot2::ylim(c(0,y_lim))
                
                plot2<-ggplot2::ggplot(plot2_data)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title2)),atop(bold(.(subtit2))))))+
                        ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                        ggplot2::ylab(label="Rates")+
                        ggplot2::xlab(label="Age (years)")+
                        ggplot2::ylim(c(0,y_lim))
                plot3<-ggplot2::ggplot(out.df)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title3)),atop(bold(.(subtit3))))))+
                        ggplot2::geom_line(ggplot2::aes(y=S,x=time))+#+theme_bw()
                        ggplot2::ylab(label="S(t)")+
                        ggplot2::xlab(label="Age (years)")#+
                #ggplot2::ylim(c(0,y_lim))
                
        } else if(constant==FALSE){
                
                #parameters=c(hb=birth_rate,hm=mortality_rate)
                if(length(birth_rate)>1 & length(mortality_rate)==1){
                        warning('Length of birth_rate is bigger than 1, assuming birth_rate is variant across time.')
                        res <- vector(length(birth_rate),mode="list")
                        for(i in seq_along(birth_rate)){
                                res[[i]]<-deSolve::ode(y=inits,times=time_var,func = equation,
                                              parms = c(hb=birth_rate[i],hm=mortality_rate))
                        }
                        plot1_data<-data.frame(age=time_var,y=birth_rate)
                        plot2_data<-data.frame(age=time_var,y=rep(mortality_rate,times=length(time_var)))
                        
                } else if(length(birth_rate)==1 & length(mortality_rate)>1){
                        warning('Length of mortality_rate is bigger than 1, assuming mortality_rate is variant across time.')
                        res <- vector(length(mortality_rate),mode="list")
                        for(i in seq_along(mortality_rate)){
                                res[[i]]<-deSolve::ode(y=inits,times=time_var,func = equation,
                                              parms = c(hb=birth_rate,#[i],
                                                        hm=mortality_rate[i]))
                        }
                        plot1_data<-data.frame(age=time_var,y=rep(birth_rate,times=length(time_var)) )
                        plot2_data<-data.frame(age=time_var,y=mortality_rate) #rep(mortality_rate,times=length(time_var)))
                } else if(length(birth_rate)>1 & length(mortality_rate)>1){
                        warning('Lengths of birth_rate and mortality_rate are bigger than 1, assuming birth_rate and mortality_rate are variant across time.')
                        res <- vector(length(mortality_rate),mode="list")
                        for(i in seq_along(mortality_rate)){
                                res[[i]]<-deSolve::ode(y=inits,times=time_var,func = equation,
                                              parms = c(hb=birth_rate[i],
                                                        hm=mortality_rate[i]))
                        }
                        plot1_data<-data.frame(age=time_var,y=birth_rate)
                        plot2_data<-data.frame(age=time_var,y=mortality_rate)
                }
                
                cout_index<-1
                new_df<-data.frame()
                for(i in 1:length(res)){
                        new_df<-rbind(new_df,
                                      res[[i]][cout_index,])
                        cout_index<-cout_index+1
                }
                
                colnames(new_df)<-c('time','S')
                
                if(max(birth_rate)>max(mortality_rate)){
                        y_lim=max(birth_rate)
                } else {
                        y_lim=max(mortality_rate)
                }
                
                title1 <- bquote("Birth rate")
                #subtit1 <- bquote(list(hb==.(parameters[1])#,~gamma==.(parameters[2])
                #))
                
                
                title2 <- bquote("Mortality rate")
                #subtit2 <- bquote(list(hm==.(parameters[2])#,~gamma==.(parameters[2])
                #))
                
                
                title3 <- bquote("Single compartmental model")
                #subtit3 <- bquote(list(hb==.(parameters[1]),~hm==.(parameters[2])
                #))
                
                #ggplot(new_df,ggplot2::aes(x=time,y=S))+
                #        ggplot2::geom_line()
                plot1<-ggplot2::ggplot(plot1_data)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title1))#,atop(bold(.(subtit1)))
                        )
                        )
                        )+
                        ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                        ggplot2::ylab(label="Rates")+
                        ggplot2::xlab(label="Age (years)")+
                        ggplot2::ylim(c(0,y_lim))
                
                plot2<-ggplot2::ggplot(plot2_data)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title2))#,atop(bold(.(subtit2)))
                        )
                        )
                        )+
                        ggplot2::geom_line(ggplot2::aes(y=y,x=age))+#+theme_bw()
                        ggplot2::ylab(label="Rates")+
                        ggplot2::xlab(label="Age (years)")+
                        ggplot2::ylim(c(0,y_lim))
                plot3<-ggplot2::ggplot(new_df)+
                        ggplot2::ggtitle(bquote(atop(bold(.(title3))#,atop(bold(.(subtit3)))
                        )
                        )
                        )+
                        ggplot2::geom_line(ggplot2::aes(y=S,x=time))+#+theme_bw()
                        ggplot2::ylab(label="S(t)")+
                        ggplot2::xlab(label="Age (years)")        
                
        } else {
                stop('constant should be defined!')
        }
        
        
        #hb_vector<-
        #hm_vector<-seq(0.01,0.015,length.out = length(time_var))
        
        #hb_res <- vector(length(hb_vector),mode="list")
        #for(i in seq_along(hb_vector)){
        #        hb_res[[i]]<-ode(y=inits,times=time_var,func = equation,
        #                         parms = c(hb=hb_vector[i],hm=hm_vector[i]))
        #}
        #names(hb_res) <- hb_vector  ## to get beta value incorporated in results
        #dd <- dplyr::bind_rows(lapply(hb_res,as.data.frame),.id="hb")
        #dd$hb <- as.numeric(dd$hb)
        
        #require(gridExtra)
        final_plots<-gridExtra::grid.arrange(plot1,plot3,plot2,layout_matrix=matrix(c(1,2,3,2),byrow = TRUE,ncol=2))
        
        return(final_plots)
}