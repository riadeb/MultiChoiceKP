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

public class Plot extends JFrame {
    public Plot(String title, LinkedList<Term> [] data, int _class){
        super(title);
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries series = new XYSeries("_class");
        for (Term pp :data[_class]) {
            series.add(pp.weight, pp.profit);

        }
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createScatterPlot(
                "Channel "+ String.valueOf(_class),
                "X-Axis",
                "Y-Axis",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);

    }
}
