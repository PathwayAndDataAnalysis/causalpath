import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.ExperimentData;

public class DummyChangeDetector implements OneDataChangeDetector {
    double changeValue;

    public DummyChangeDetector(double changeValue){
        this.changeValue = changeValue;
    }
    @Override
    public int getChangeSign(ExperimentData data) {
        return 0;
    }

    @Override
    public double getChangeValue(ExperimentData data) {
        return changeValue;
    }

    @Override
    public OneDataChangeDetector makeACopy() {
        return null;
    }

}
