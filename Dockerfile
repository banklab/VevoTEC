FROM r-base:4.3.3

RUN R -e "install.packages(c('shiny', 'plotly'))"

CMD R main.R
