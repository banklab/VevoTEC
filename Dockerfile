FROM rocker/shiny-verse:4.0.3

RUN R -e "install.packages('plotly')"

CMD R main.R
