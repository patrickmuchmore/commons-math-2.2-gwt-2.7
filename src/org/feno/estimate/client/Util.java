package org.feno.estimate.client;

import com.google.gwt.dom.client.Style.Cursor;
import com.google.gwt.user.client.ui.RootPanel;

public class Util {

	public native static final void debugger() /*-{
												debugger;
												}-*/;

	public static final void showWaitCursor() {
		RootPanel.getBodyElement().getStyle().setCursor(Cursor.WAIT);
	}

	public static final void showDefaultCursor() {
		RootPanel.getBodyElement().getStyle().setCursor(Cursor.DEFAULT);
	}

	public static final String getMessage(Throwable throwable) {
		String ret = "";
		while (throwable != null) {
			if (throwable instanceof com.google.gwt.event.shared.UmbrellaException) {
				for (Throwable thr2 : ((com.google.gwt.event.shared.UmbrellaException) throwable).getCauses()) {
					if (ret != "")
						ret += "\nCaused by: ";
					ret += thr2.toString();
					ret += "\n  at " + getMessage(thr2);
				}
			} else if (throwable instanceof com.google.web.bindery.event.shared.UmbrellaException) {
				for (Throwable thr2 : ((com.google.web.bindery.event.shared.UmbrellaException) throwable).getCauses()) {
					if (ret != "")
						ret += "\nCaused by: ";
					ret += thr2.toString();
					ret += "\n  at " + getMessage(thr2);
				}
			} else {
				if (ret != "")
					ret += "\nCaused by: ";
				ret += throwable.toString();
				for (StackTraceElement sTE : throwable.getStackTrace())
					ret += "\n  at " + sTE;
			}
			throwable = throwable.getCause();
		}
		return ret;
	}

}
