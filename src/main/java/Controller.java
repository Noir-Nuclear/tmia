import javafx.collections.FXCollections;
import javafx.fxml.Initializable;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;

import java.net.URL;
import java.util.List;
import java.util.ResourceBundle;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

public class Controller implements Initializable {
    public ScatterChart pointsChart;
    public ComboBox matricesDictionary;
    public TextArea matrixField;
    public TextField kField;
    public Button rerenderChart;
    public TextArea resultField;
    /**
     * Матрицы преобразований
     */
    private List<List<List<Double>>> matrices;

    private double[] interval = {0, 2 * Math.PI};

    @Override
    public void initialize(URL location, ResourceBundle resources) {
        initializeMatrices();
        matricesDictionary.setItems(FXCollections.observableArrayList(List.of(1, 2, 3, 4, 5, 6, 7, 8)));
        matricesDictionary.setValue(1);
        fillMatrixField();
        matricesDictionary.setOnAction(event -> fillMatrixField());
        List<List<Double>> points = DoubleStream.iterate(
                interval[0], t -> t < interval[1], t -> t + 0.01).
                mapToObj(t -> List.of(fX(t), fY(t))).
                collect(Collectors.toList());
        renderChart(points, "исходное");
        rerenderChart.setOnAction(event -> solve());
    }

    private void solve() {
        int k;
        try {
            k = Integer.parseInt(kField.getText());
        } catch (Exception e) {
            resultField.setText("Заполните поле масштаба!");
            return;
        }
        pointsChart.getData().clear();
        int indexOfMatrix = (int) matricesDictionary.getValue() - 1;
        drawNominalLine(indexOfMatrix, k);
        double beta = calculateBeta(indexOfMatrix);
        double alpha = calculateAlpha(indexOfMatrix);
        double alpha1 = calculateAlpha1(beta, indexOfMatrix);
        double alpha2 = calculateAlpha2(beta, indexOfMatrix);

        List<Double> kx = calculateKx(alpha, beta, indexOfMatrix);
        List<Double> ky = calculateKy(alpha, beta, indexOfMatrix);

        renderChart(
                performToOrtoBasis(alpha, beta, kx.get(0), ky.get(0), k),
                "итоговое"
        );
        resultField.setText(
                "kx1:   " + kx.get(0) + ", kx2: " + kx.get(1) + "\n" +
                        "ky1:   " + ky.get(0) + ", ky2: " + ky.get(1) + "\n" +
                        "alpha:    " + alpha * 180 / Math.PI + "\n" +
                        "alpha1:    " + alpha1 * 180 / Math.PI + "\n" +
                        "alpha2:    " + alpha2 * 180 / Math.PI + "\n" +
                        "beta:  " + beta * 180 / Math.PI
        );
    }

    private List<List<Double>> performToOrtoBasis(double alpha, double beta, double kx, double ky, int k) {
        List<List<Double>> points = DoubleStream.iterate(
                interval[0], t -> t < interval[1], t -> t + 0.01).
                mapToObj(t -> List.of(fX(t), fY(t))).
                collect(Collectors.toList());
        double A = 1.0 / Math.sqrt(
                Math.pow(Math.cos(alpha), 2.0) +
                        Math.pow(Math.sin(-beta), 2.0)
        );
        double B = 1.0 / Math.sqrt(
                Math.pow(Math.sin(alpha), 2) +
                        Math.pow(Math.cos(beta), 2)
        );
        points = points.stream().map(point -> List.of(
                A * point.get(0),
                B * point.get(1)
        )).collect(Collectors.toList());
        if (false) {
            points = points.stream().map(point -> List.of(
                    (point.get(0) * Math.cos(alpha) +
                            point.get(1) * Math.sin(alpha)) *
                            kx * k,
                    (point.get(0) * Math.cos(alpha) +
                            point.get(1) * Math.sin(alpha)) *
                            ky * k
            )).collect(Collectors.toList());
        } else {
            points = points.stream().map(point -> List.of(
                    (point.get(0) * -Math.cos(alpha) +
                            point.get(1) * Math.sin(alpha)) *
                            kx * k,
                    (point.get(0) * Math.sin(beta) +
                            point.get(1) * Math.cos(beta)) *
                            ky * k
            )).collect(Collectors.toList());
        }
        points = points.stream().map(point -> List.of(
                point.get(0) * Math.cos(alpha) - point.get(1) * Math.sin(alpha),
                point.get(0) * Math.sin(alpha) + point.get(1) * Math.cos(alpha)
        )).collect(Collectors.toList());
        return points;
    }

    private void fillMatrixField() {
        List<List<Double>> matrix = matrices.get((Integer) matricesDictionary.getValue() - 1);
        matrixField.setText(
                matrix.get(0).get(0) + "  " + matrix.get(0).get(1) + "\n" +
                        matrix.get(1).get(0) + "    " + matrix.get(1).get(1)
        );
    }

    private double calculateBeta(int i) {
        double a = matrices.get(i).get(0).get(0);
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        double b = matrices.get(i).get(1).get(1);
        return Math.atan2(
                2.0 * (a * h + b * g),
                ((Math.pow(a, 2.0) - Math.pow(h, 2.0)) - (Math.pow(b, 2.0) - Math.pow(g, 2.0)))
        ) / 2;
    }

    private double calculateAlpha(int i) {
        double a = matrices.get(i).get(0).get(0);
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        double b = matrices.get(i).get(1).get(1);
        return Math.atan2(
                2.0 * (b * h + a * g),
                ((Math.pow(a, 2.0) + Math.pow(h, 2.0)) - (Math.pow(b, 2.0) + Math.pow(g, 2.0)))
        ) / 2;
    }

    private double calculateAlpha1(double beta, int i) {
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        return Math.atan2(
                Math.sin(beta) - h * Math.cos(beta),
                Math.cos(beta) - g * Math.sin(beta)
        );
    }

    private double calculateAlpha2(double beta, int i) {
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        return Math.atan2(
                Math.sin(beta) + g * Math.cos(beta),
                Math.cos(beta) + h * Math.sin(beta)
        );
    }

    private List<Double> calculateKx(double alpha, double beta, int i) {
        double a = matrices.get(i).get(0).get(0);
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        double b = matrices.get(i).get(1).get(1);
        return List.of(
                (a * Math.cos(beta) + h * Math.sin(beta)) / Math.cos(alpha),
                (b * Math.sin(beta) + g * Math.cos(beta)) / Math.sin(alpha)
        );
    }

    private List<Double> calculateKy(double alpha, double beta, int i) {
        double a = matrices.get(i).get(0).get(0);
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        double b = matrices.get(i).get(1).get(1);
        return List.of(
                (a * Math.sin(beta) - h * Math.cos(beta)) / Math.sin(alpha),
                (b * Math.cos(beta) - g * Math.sin(beta)) / Math.cos(alpha)
        );
    }

    private void drawNominalLine(int i, int k) {
        double a = matrices.get(i).get(0).get(0);
        double h = matrices.get(i).get(0).get(1);
        double g = matrices.get(i).get(1).get(0);
        double b = matrices.get(i).get(1).get(1);
        List<List<Double>> points = DoubleStream.iterate(
                interval[0], t -> t < interval[1], t -> t + 0.01).
                mapToObj(t -> List.of(fXS(t, a, h, k), fYS(t, g, b, k))).
                collect(Collectors.toList());
        renderChart(points, "номинальное");
    }

    private void initializeMatrices() {
        matrices = List.of(
                List.of(
                        List.of(0.01, -1.0),
                        List.of(1.0, 0.7)
                ),//1
                List.of(
                        List.of(0.01, 1.0),
                        List.of(-1.0, 0.7)
                ),//2
                List.of(
                        List.of(0.01, -1.0),
                        List.of(-1.0, 0.7)
                ),//3
                List.of(
                        List.of(0.01, 1.0),
                        List.of(1.0, 0.7)
                ),//4
                List.of(
                        List.of(-0.01, -1.0),
                        List.of(1.0, 0.7)
                ),//5
                List.of(
                        List.of(-0.01, 1.0),
                        List.of(-1.0, 0.7)
                ),//6
                List.of(
                        List.of(-0.01, -1.0),
                        List.of(-1.0, 0.7)
                ),//7
                List.of(
                        List.of(-0.01, 1.0),
                        List.of(1.0, 0.7)
                )//8
        );

    }

    private double fX(double t) {
        return (1 + 2 * Math.cos(t)) * Math.sin(t);
    }

    private double fY(double t) {
        return (1 + 2 * Math.cos(t)) * Math.cos(t);
    }

    private double fXS(double t, double a, double h, double k) {
        return k * (a * fX(t) + h * fY(t));
    }

    private double fYS(double t, double g, double b, double k) {
        return k * (g * fX(t) + b * fY(t));
    }

    private void renderChart(List<List<Double>> points, String legend) {
        XYChart.Series series = new XYChart.Series();
        series.setName(legend);
        points.forEach(point ->
                series.getData().add(new XYChart.Data(point.get(0), point.get(1))));
        pointsChart.getData().add(series);
    }
}
