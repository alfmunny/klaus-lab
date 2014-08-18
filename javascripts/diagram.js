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
              text: 'Training Frequency'
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
            //layout: 'horizontal',
            layout: 'vertical',
            align: 'left',
            x: 60,
            y: -60,
            //align: '',
            //verticalAlign: 'middle',
            borderWidth: 2,
            floating: true
        },
        series: data
    });
});
