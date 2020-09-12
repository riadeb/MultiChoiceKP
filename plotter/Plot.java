package plotter;

import containers.Term;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import java.util.LinkedList;

import javax.swing.*;
import java.util.ArrayList;

public class Plot extends JFrame {
    public Plot(String title, LinkedList<Term> [] data, int channel){
        super(title);
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries series = new XYSeries("channel");
        for (Term pp :data[channel]) {
            series.add(pp.p, pp.r);

        }
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createScatterPlot(
                "Channel "+ String.valueOf(channel),
                "X-Axis",
                "Y-Axis",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);

    }
}
