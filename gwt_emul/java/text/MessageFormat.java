package java.text;

import java.util.Locale;

public class MessageFormat extends Format {

  public MessageFormat(String localizedString, Locale locale) {
    // TODO Auto-generated constructor stub
  }

  @Override
  public StringBuffer format(Object object, StringBuffer buffer,
      FieldPosition field) {
    // FIXME 
    return buffer.append(object.toString());
  }

  @Override
  public Object parseObject(String string, ParsePosition position) {
    // TODO Auto-generated method stub
    return null;
  }

}
