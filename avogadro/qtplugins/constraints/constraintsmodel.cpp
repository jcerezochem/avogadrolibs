/****************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/

#include "constraintsmodel.h"

#include <QtCore/QDebug>
#include <QtCore/QString>

using namespace std;

namespace Avogadro {
using Core::Constraint;

namespace QtPlugins {

void ConstraintsModel::emitDataChanged()
{
  emit dataChanged(QModelIndex(), QModelIndex());
}

void ConstraintsModel::setConstraints(
  const std::vector<Core::Constraint>& constraints)
{
  beginResetModel();
  m_constraints = constraints;
  endResetModel();
}

int ConstraintsModel::rowCount(const QModelIndex&) const
{
  return m_constraints.size();
}

int ConstraintsModel::columnCount(const QModelIndex&) const
{
  // Type, value, atom 1, 2, 3, 4, scan initial, scan end, scan steps
  return 9;
}

QVariant ConstraintsModel::formattedConstraintValue(
  const Core::Constraint& constraint, int role) const
{
  if (role == Qt::UserRole || role == Qt::EditRole)
    return constraint.value();

  if (constraint.type() == 1)
    return tr("%L1 Å").arg(constraint.value(), 0, 'f', 3);
  if (constraint.type() == 2 || constraint.type() == 3)
    return tr("%L1 °").arg(constraint.value(), 0, 'f', 3);

  return "--";
}

QVariant ConstraintsModel::data(const QModelIndex& index, int role) const
{
  if (!index.isValid())
    return QVariant();

  if (index.row() >= static_cast<int>(m_constraints.size()))
    return QVariant();

  const Constraint& currentConstraint = m_constraints[index.row()];
  Index aIndex = currentConstraint.aIndex();
  Index bIndex = currentConstraint.bIndex();
  Index cIndex = currentConstraint.cIndex();
  Index dIndex = currentConstraint.dIndex();

  if (role == Qt::DisplayRole || role == Qt::UserRole || role == Qt::EditRole)
    switch (index.column()) {
      case 0:
        if (currentConstraint.type() == 1)
          return tr("Distance");
        else if (currentConstraint.type() == 2)
          return tr("Angle");
        else if (currentConstraint.type() == 3)
          return tr("Torsion Angle");
        else if (currentConstraint.type() == 4)
          return tr("Freeze Atom", "fix / remain constant");
        else if (currentConstraint.type() == 5)
          return tr("Freeze X Axis", "fix / remain constant");
        else if (currentConstraint.type() == 6)
          return tr("Freeze Y Axis", "fix / remain constant");
        else if (currentConstraint.type() == 7)
          return tr("Freeze Z Axis", "fix / remain constant");
        break;
      case 1:
        return formattedConstraintValue(currentConstraint, role);
      case 2:
        if (aIndex != MaxIndex)
          return QVariant(static_cast<qulonglong>(aIndex) + 1);
        else
          return "--";
      case 3:
        if (bIndex != MaxIndex)
          return QVariant(static_cast<qulonglong>(bIndex) + 1);
        else
          return "--";
      case 4:
        if (cIndex != MaxIndex)
          return QVariant(static_cast<qulonglong>(cIndex) + 1);
        else
          return "--";
      case 5:
        if (dIndex != MaxIndex)
          return QVariant(static_cast<qulonglong>(dIndex) + 1);
        else
          return "--";
      case 6:
        if (role == Qt::EditRole)
          return currentConstraint.hasScan() ? currentConstraint.scanInitial()
                                             : currentConstraint.value();
        if (role == Qt::UserRole)
          return currentConstraint.hasScan() ? currentConstraint.scanInitial()
                                             : QVariant();
        if (currentConstraint.hasScan())
          return formattedConstraintValue(
            Constraint(aIndex, bIndex, cIndex, dIndex,
                       currentConstraint.scanInitial()),
            Qt::DisplayRole);
        return "--";
      case 7:
        if (role == Qt::EditRole)
          return currentConstraint.hasScan() ? currentConstraint.scanEnd()
                                             : currentConstraint.value();
        if (role == Qt::UserRole)
          return currentConstraint.hasScan() ? currentConstraint.scanEnd()
                                             : QVariant();
        if (currentConstraint.hasScan())
          return formattedConstraintValue(
            Constraint(aIndex, bIndex, cIndex, dIndex, currentConstraint.scanEnd()),
            Qt::DisplayRole);
        return "--";
      case 8:
        if (role == Qt::EditRole || role == Qt::UserRole)
          return currentConstraint.scanSteps();
        if (currentConstraint.hasScan())
          return currentConstraint.scanSteps();
        return "--";
    }

  return QVariant();
}

bool ConstraintsModel::updateScanData(Core::Constraint& constraint, int column,
                                      const QVariant& value)
{
  bool ok = false;
  if (column == 8) {
    int steps = value.toInt(&ok);
    if (!ok)
      return false;

    if (steps < 2) {
      constraint.clearScan();
      return true;
    }

    Real initial = constraint.hasScan() ? constraint.scanInitial()
                                        : constraint.value();
    Real end = constraint.hasScan() ? constraint.scanEnd() : constraint.value();
    constraint.setScan(initial, end, steps);
    return true;
  }

  Real editedValue = value.toDouble(&ok);
  if (!ok)
    return false;

  Real initial = constraint.hasScan() ? constraint.scanInitial()
                                      : constraint.value();
  Real end = constraint.hasScan() ? constraint.scanEnd() : constraint.value();
  int steps = constraint.hasScan() ? constraint.scanSteps() : 2;

  if (column == 6)
    initial = editedValue;
  else if (column == 7)
    end = editedValue;
  else
    return false;

  constraint.setScan(initial, end, steps);
  return true;
}

bool ConstraintsModel::setData(const QModelIndex& index, const QVariant& value,
                               int role)
{
  if (role != Qt::EditRole || !index.isValid() ||
      index.row() >= static_cast<int>(m_constraints.size())) {
    return false;
  }

  if (index.column() < 6 || index.column() > 8)
    return false;

  auto& constraint = m_constraints[index.row()];
  if (!updateScanData(constraint, index.column(), value))
    return false;

  emit dataChanged(index, index, { Qt::DisplayRole, Qt::EditRole, Qt::UserRole });
  return true;
}

Qt::ItemFlags ConstraintsModel::flags(const QModelIndex& index) const
{
  if (!index.isValid())
    return Qt::NoItemFlags;

  Qt::ItemFlags itemFlags = QAbstractTableModel::flags(index);
  if (index.column() >= 6 && index.column() <= 8)
    itemFlags |= Qt::ItemIsEditable;

  return itemFlags;
}

QVariant ConstraintsModel::headerData(int section, Qt::Orientation orientation,
                                      int role) const
{
  if (role != Qt::DisplayRole)
    return QVariant();

  if (orientation == Qt::Horizontal) {
    switch (section) {
      case 0:
        return tr("Type");
      case 1:
        return tr("Value");
      case 2:
        return tr("Atom 1");
      case 3:
        return tr("Atom 2");
      case 4:
        return tr("Atom 3");
      case 5:
        return tr("Atom 4");
      case 6:
        return tr("Start");
      case 7:
        return tr("End");
      case 8:
        return tr("Steps");
    }
  }

  return section + 1;
}

void ConstraintsModel::addConstraint(int type, int a, int b, int c, int d,
                                     double value)
{
  beginInsertRows(QModelIndex(), m_constraints.size(), m_constraints.size());
  m_constraints.push_back(Constraint(a, b, c, d, value));
  m_constraints.back().setType(static_cast<Constraint::Type>(type));
  endInsertRows();
}


void ConstraintsModel::setConstraint(int index,
                                     const Core::Constraint& constraint)
{
  if (index < 0 || index >= static_cast<int>(m_constraints.size()))
    return;

  m_constraints[index] = constraint;
  emit dataChanged(this->index(index, 0), this->index(index, columnCount() - 1),
                   { Qt::DisplayRole, Qt::EditRole, Qt::UserRole });
}

void ConstraintsModel::clear()
{
  if (m_constraints.size()) {
    beginRemoveRows(QModelIndex(), 0, m_constraints.size() - 1);
    m_constraints.clear();
    endRemoveRows();
  }
}

void ConstraintsModel::deleteConstraint(int index)
{
  if (m_constraints.size() && (index >= 0)) {
    beginRemoveRows(QModelIndex(), index, index);
    auto position = m_constraints.begin() + index;
    m_constraints.erase(position);
    endRemoveRows();
  }
}

Core::Constraint ConstraintsModel::constraint(int index)
{
  if (index < 0 || index >= static_cast<int>(m_constraints.size()))
    return Constraint(MaxIndex, MaxIndex, MaxIndex, MaxIndex, 0.0);
  else
    return m_constraints[index];
}

} // namespace QtPlugins
} // end namespace Avogadro
