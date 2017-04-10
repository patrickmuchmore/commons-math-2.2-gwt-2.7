/* 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package java.text;

import java.io.InvalidObjectException;
import java.util.Locale;

public abstract class NumberFormat extends Format {

  private static final long serialVersionUID = -2308460125733713944L;

  /**
   * Field constant identifying the integer part of a number.
   */
  public static final int INTEGER_FIELD = 0;

  /**
   * Field constant identifying the fractional part of a number.
   */
  public static final int FRACTION_FIELD = 1;

  private boolean groupingUsed = true, parseIntegerOnly = false;

  private int maximumIntegerDigits = 40, minimumIntegerDigits = 1,
      maximumFractionDigits = 3, minimumFractionDigits = 0;

  /**
   * Constructs a new instance of {@code NumberFormat}.
   */
  public NumberFormat() {
  }

  /**
   * Indicates whether this number format formats and parses numbers using a
   * grouping separator.
   * 
   * @return {@code true} if a grouping separator is used; {@code false}
   *         otherwise.
   */
  public boolean isGroupingUsed() {
    return groupingUsed;
  }

  /**
   * Sets whether this number format formats and parses numbers using a
   * grouping separator.
   * 
   * @param value
   *            {@code true} if a grouping separator is used; {@code false}
   *            otherwise.
   */
  public void setGroupingUsed(boolean value) {
    groupingUsed = value;
  }

  @Override
  public StringBuffer format(Object object, StringBuffer buffer,
      FieldPosition field) {
    // FIXME 
    return buffer.append(object.toString());
  }

  public static Locale[] getAvailableLocales() {
    // FIXME
    return new Locale[]{Locale.getDefault()};
  }

  /**
   * Parses a {@code Number} from the specified string using the rules of this
   * number format.
   * 
   * @param string
   *            the string to parse.
   * @return the {@code Number} resulting from the parsing.
   * @throws ParseException
   *            if an error occurs during parsing.
   */
  public Number parse(String string) throws ParseException {
    ParsePosition pos = new ParsePosition(0);
    Number number = parse(string, pos);
    if (pos.getIndex() == 0) {
      // text.1D=Unparseable number: {0}
      throw new ParseException(
          Messages.getString("text.1D", string), pos.getErrorIndex()); //$NON-NLS-1$
    }
    return number;
  }

  /**
   * Parses a {@code Number} from the specified string starting at the index
   * specified by {@code position}. If the string is successfully parsed then
   * the index of the {@code ParsePosition} is updated to the index following
   * the parsed text. On error, the index is unchanged and the error index of
   * {@code ParsePosition} is set to the index where the error occurred.
   * 
   * @param string
   *            the string to parse.
   * @param position
   *            input/output parameter, specifies the start index in
   *            {@code string} from where to start parsing. If parsing is
   *            successful, it is updated with the index following the parsed
   *            text; on error, the index is unchanged and the error index is
   *            set to the index where the error occurred.
   * @return the {@code Number} resulting from the parse or {@code null} if
   *         there is an error.
   */
  public abstract Number parse(String string, ParsePosition position);

  @Override
  public final Object parseObject(String string, ParsePosition position) {
    if (position == null) {
      // text.1A=position is null
      throw new NullPointerException(Messages.getString("text.1A")); //$NON-NLS-1$
    }

    try {
      return parse(string, position);
    } catch (Exception e) {
      return null;
    }
  }

  public abstract StringBuffer format(long number, StringBuffer buffer, FieldPosition position);

  public abstract StringBuffer format(double number, StringBuffer buffer, FieldPosition position);

  public final static NumberFormat getInstance(Locale locale) {
    return getNumberInstance(locale);
  }

  public final static NumberFormat getNumberInstance(Locale locale) {
    // TODO Auto-generated method stub
    return null;
  }
  
  public final static NumberFormat getPercentInstance() {
    return getPercentInstance(Locale.getDefault());
  }

  public final static NumberFormat getPercentInstance(Locale locale) {
    // TODO Auto-generated method stub
    return null;
  }
  
  /**
   * Sets the maximum number of fraction digits that are printed when
   * formatting. If the maximum is less than the number of fraction digits,
   * the least significant digits are truncated.
   * 
   * @param value
   *            the maximum number of fraction digits.
   */
  public void setMaximumFractionDigits(int value) {
    maximumFractionDigits = value < 0 ? 0 : value;
    if (maximumFractionDigits < minimumFractionDigits) {
      minimumFractionDigits = maximumFractionDigits;
    }
  }

  /**
   * Sets the new maximum count of integer digits that are printed when
   * formatting. If the maximum is less than the number of integer digits, the
   * most significant digits are truncated.
   * 
   * @param value
   *            the new maximum number of integer numerals for display.
   */
  public void setMaximumIntegerDigits(int value) {
    maximumIntegerDigits = value < 0 ? 0 : value;
    if (maximumIntegerDigits < minimumIntegerDigits) {
      minimumIntegerDigits = maximumIntegerDigits;
    }
  }

  /**
   * Specifies if this number format should parse numbers only as integers or
   * else as any kind of number. If this method is called with a {@code true}
   * value then subsequent parsing attempts will stop if a decimal separator
   * is encountered.
   * 
   * @param value
   *            {@code true} to only parse integers, {@code false} to parse
   *            integers as well as fractions.
   */
  public void setParseIntegerOnly(boolean value) {
    parseIntegerOnly = value;
  }

  /**
   * The instances of this inner class are used as attribute keys and values
   * in {@code AttributedCharacterIterator} that the
   * {@link NumberFormat#formatToCharacterIterator(Object)} method returns.
   * <p>
   * There is no public constructor in this class, the only instances are the
   * constants defined here.
   * <p>
   */
  public static class Field extends Format.Field {

    private static final long serialVersionUID = 7494728892700160890L;

    /**
     * This constant stands for the number sign.
     */
    public static final Field SIGN = new Field("sign"); //$NON-NLS-1$

    /**
     * This constant stands for the integer part of the number.
     */
    public static final Field INTEGER = new Field("integer"); //$NON-NLS-1$

    /**
     * This constant stands for the fraction part of the number.
     */
    public static final Field FRACTION = new Field("fraction"); //$NON-NLS-1$

    /**
     * This constant stands for the exponent part of the number.
     */
    public static final Field EXPONENT = new Field("exponent"); //$NON-NLS-1$

    /**
     * This constant stands for the exponent sign symbol.
     */
    public static final Field EXPONENT_SIGN = new Field("exponent sign"); //$NON-NLS-1$

    /**
     * This constant stands for the exponent symbol.
     */
    public static final Field EXPONENT_SYMBOL = new Field("exponent symbol"); //$NON-NLS-1$

    /**
     * This constant stands for the decimal separator.
     */
    public static final Field DECIMAL_SEPARATOR = new Field(
        "decimal separator"); //$NON-NLS-1$

    /**
     * This constant stands for the grouping separator.
     */
    public static final Field GROUPING_SEPARATOR = new Field(
        "grouping separator"); //$NON-NLS-1$

    /**
     * This constant stands for the percent symbol.
     */
    public static final Field PERCENT = new Field("percent"); //$NON-NLS-1$

    /**
     * This constant stands for the permille symbol.
     */
    public static final Field PERMILLE = new Field("per mille"); //$NON-NLS-1$

    /**
     * This constant stands for the currency symbol.
     */
    public static final Field CURRENCY = new Field("currency"); //$NON-NLS-1$

    /**
     * Constructs a new instance of {@code NumberFormat.Field} with the
     * given field name.
     *
     * @param fieldName
     *            the field name.
     */
    protected Field(String fieldName) {
      super(fieldName);
    }

    /**
     * Resolves instances that are deserialized to the constant
     * {@code NumberFormat.Field} values.
     *
     * @return the resolved field object.
     * @throws InvalidObjectException
     *             if an error occurs while resolving the field object.
     */
    @Override
    protected Object readResolve() throws InvalidObjectException {
      if (this.equals(INTEGER)) {
        return INTEGER;
      }
      if (this.equals(FRACTION)) {
        return FRACTION;
      }
      if (this.equals(EXPONENT)) {
        return EXPONENT;
      }
      if (this.equals(EXPONENT_SIGN)) {
        return EXPONENT_SIGN;
      }
      if (this.equals(EXPONENT_SYMBOL)) {
        return EXPONENT_SYMBOL;
      }
      if (this.equals(CURRENCY)) {
        return CURRENCY;
      }
      if (this.equals(DECIMAL_SEPARATOR)) {
        return DECIMAL_SEPARATOR;
      }
      if (this.equals(GROUPING_SEPARATOR)) {
        return GROUPING_SEPARATOR;
      }
      if (this.equals(PERCENT)) {
        return PERCENT;
      }
      if (this.equals(PERMILLE)) {
        return PERMILLE;
      }
      if (this.equals(SIGN)) {
        return SIGN;
      }
      // text.02=Unknown attribute
          throw new InvalidObjectException(Messages.getString("text.02")); //$NON-NLS-1$
    }
  }

  public NumberFormat clone() {
    // TODO Auto-generated method stub
    return null;
  }

}
