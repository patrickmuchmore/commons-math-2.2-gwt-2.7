package org.feno.estimate.client.mcmc;

import java.util.Collection;
import java.util.HashMap;

import org.apache.commons.math.random.RandomGenerator;
import org.feno.estimate.client.model.Model;
import org.feno.estimate.client.model.Parameter;

import com.google.gwt.typedarrays.client.Float64ArrayNative;

public class SamplePosterior {

	public static HashMap<String, Float64ArrayNative> sample(Model model, int steps, RandomGenerator rng) {

		Collection<Parameter> parameters = model.getParameters();
		HashMap<String, Float64ArrayNative> parameterSamples = new HashMap<String, Float64ArrayNative>(
				parameters.size());

		for (Parameter param : parameters) {
			parameterSamples.put(param.getName(), Float64ArrayNative.create(steps));
		}

		HashMap<String, Double> currentParameterValues = new HashMap<String, Double>(parameters.size());
		double currentLogLikelihood = model.logLikelihood();
		double proposedLogLikelihood = Double.NEGATIVE_INFINITY;

		for (int i = 0; i < steps; i++) {

			for (Parameter param : parameters) {
				currentParameterValues.put(param.getName(), param.getValue());
				param.setValue(param.reflectProposal(param.getValue() + rng.nextGaussian() * param.getPropSd()));
			}

			proposedLogLikelihood = model.logLikelihood();

			if (Math.log(rng.nextDouble()) < proposedLogLikelihood - currentLogLikelihood) {
				currentLogLikelihood = proposedLogLikelihood;
			} else {
				for (Parameter param : parameters) {
					param.setValue(currentParameterValues.get(param.getName()));
				}
			}

			for (Parameter param : parameters) {
				parameterSamples.get(param.getName()).set(i, param.getValue());
			}
		}

		return parameterSamples;
	}

}
