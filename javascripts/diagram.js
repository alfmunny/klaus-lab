$(function () {
    $('#diagram').highcharts({
        title: {
            text: 'Recognition Acurracy',
            x: -20 //center
        },
        subtitle: {
            text: 'For all models',
            x: -20
        },
        xAxis: {
            title: {
              text: 'Training Frequecy'
            },

            categories: ['1', '2', '3', '4', '5', '6',
                '7', '8', '9', '10', '11', '12',
                '13', '14', '15', '16', '17', '18', '19', '20']
        },
        yAxis: {
            title: {
                text: 'Acurracy (%)'
            },
            plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
            }]
        },
        tooltip: {
            valueSuffix: '%'
        },
        legend: {
            layout: 'vertical',
            align: '',
            verticalAlign: 'middle',
            borderWidth: 0
        },
        series: [{
            name: '10 dim slant bakis',
            data: [42.81, 42.7, 42.48, 42.44, 42.67, 43.00, 42.85, 43.26, 43.52, 43.15, 26.44, 43.48,
            33.56, 37.59, 39.56, 41.11, 41.89, 41.96, 42.48, 42.74]
        }, {
            name: '10 dim slant nobakis',
            data: [38.63, 38.96, 39.37, 39.56, 40.26, 40.48, 39.74, 40.26, 40.37, 40.74, 26.15, 40.41,
            32.26, 35.04, 37.44, 38.63, 38.22, 38.89, 38.78, 38.89]
        }, {
            name: '10 dim noslant bakis',
            data: [44.37, 45.00, 45.11, 45.41, 44.70, 44.56, 44.81, 44.96, 45.26, 45.04, 36.41, 45.26,
            41.04, 43.04, 42.56, 43.15, 43.30, 43.78, 43.41, 43.81]
        },{
            name: '10 dim noslant nobakis',
            data: [38.78, 38.70, 39.04, 39.22, 39.04, 38.93, 39.04, 39.00, 39.07, 39.52, 35.41, 39.41,
            36.48, 38.11, 38.74, 38.85, 38.44, 38.59, 39.00, 38.85]
        }]
    });
});
