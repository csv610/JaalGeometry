/****************************************************************************
** Meta object code from reading C++ file 'ShapeOpDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ShapeOpDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ShapeOpDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JShapeOpDialog_t {
    QByteArrayData data[19];
    char stringdata0[424];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JShapeOpDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JShapeOpDialog_t qt_meta_stringdata_JShapeOpDialog = {
    {
        QT_MOC_LITERAL(0, 0, 14), // "JShapeOpDialog"
        QT_MOC_LITERAL(1, 15, 20), // "applyAreaConstraints"
        QT_MOC_LITERAL(2, 36, 0), // ""
        QT_MOC_LITERAL(3, 37, 26), // "applyCocircularConstraints"
        QT_MOC_LITERAL(4, 64, 24), // "applyCoplanarConstraints"
        QT_MOC_LITERAL(5, 89, 26), // "applyCrossFieldConstraints"
        QT_MOC_LITERAL(6, 116, 26), // "applyEdgeStrainConstraints"
        QT_MOC_LITERAL(7, 143, 27), // "applyFixedLengthConstraints"
        QT_MOC_LITERAL(8, 171, 23), // "applyLaplaceConstraints"
        QT_MOC_LITERAL(9, 195, 29), // "applyParallelogramConstraints"
        QT_MOC_LITERAL(10, 225, 30), // "applyTriangleStrainConstraints"
        QT_MOC_LITERAL(11, 256, 25), // "applyRectangleConstraints"
        QT_MOC_LITERAL(12, 282, 25), // "applyClosenessConstraints"
        QT_MOC_LITERAL(13, 308, 29), // "applyFixedBoundaryConstraints"
        QT_MOC_LITERAL(14, 338, 26), // "applySimilarityConstraints"
        QT_MOC_LITERAL(15, 365, 21), // "applyRigidConstraints"
        QT_MOC_LITERAL(16, 387, 18), // "displayConstraints"
        QT_MOC_LITERAL(17, 406, 5), // "solve"
        QT_MOC_LITERAL(18, 412, 11) // "closeDialog"

    },
    "JShapeOpDialog\0applyAreaConstraints\0"
    "\0applyCocircularConstraints\0"
    "applyCoplanarConstraints\0"
    "applyCrossFieldConstraints\0"
    "applyEdgeStrainConstraints\0"
    "applyFixedLengthConstraints\0"
    "applyLaplaceConstraints\0"
    "applyParallelogramConstraints\0"
    "applyTriangleStrainConstraints\0"
    "applyRectangleConstraints\0"
    "applyClosenessConstraints\0"
    "applyFixedBoundaryConstraints\0"
    "applySimilarityConstraints\0"
    "applyRigidConstraints\0displayConstraints\0"
    "solve\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JShapeOpDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    17,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   99,    2, 0x08 /* Private */,
    3,    0,  100,    2, 0x08 /* Private */,
    4,    0,  101,    2, 0x08 /* Private */,
    5,    0,  102,    2, 0x08 /* Private */,
    6,    0,  103,    2, 0x08 /* Private */,
    7,    0,  104,    2, 0x08 /* Private */,
    8,    0,  105,    2, 0x08 /* Private */,
    9,    0,  106,    2, 0x08 /* Private */,
    10,    0,  107,    2, 0x08 /* Private */,
    11,    0,  108,    2, 0x08 /* Private */,
    12,    0,  109,    2, 0x08 /* Private */,
    13,    0,  110,    2, 0x08 /* Private */,
    14,    0,  111,    2, 0x08 /* Private */,
    15,    0,  112,    2, 0x08 /* Private */,
    16,    0,  113,    2, 0x08 /* Private */,
    17,    0,  114,    2, 0x08 /* Private */,
    18,    0,  115,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JShapeOpDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JShapeOpDialog *_t = static_cast<JShapeOpDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->applyAreaConstraints();
            break;
        case 1:
            _t->applyCocircularConstraints();
            break;
        case 2:
            _t->applyCoplanarConstraints();
            break;
        case 3:
            _t->applyCrossFieldConstraints();
            break;
        case 4:
            _t->applyEdgeStrainConstraints();
            break;
        case 5:
            _t->applyFixedLengthConstraints();
            break;
        case 6:
            _t->applyLaplaceConstraints();
            break;
        case 7:
            _t->applyParallelogramConstraints();
            break;
        case 8:
            _t->applyTriangleStrainConstraints();
            break;
        case 9:
            _t->applyRectangleConstraints();
            break;
        case 10:
            _t->applyClosenessConstraints();
            break;
        case 11:
            _t->applyFixedBoundaryConstraints();
            break;
        case 12:
            _t->applySimilarityConstraints();
            break;
        case 13:
            _t->applyRigidConstraints();
            break;
        case 14:
            _t->displayConstraints();
            break;
        case 15:
            _t->solve();
            break;
        case 16:
            _t->closeDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JShapeOpDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JShapeOpDialog.data,
        qt_meta_data_JShapeOpDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JShapeOpDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JShapeOpDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JShapeOpDialog.stringdata0))
        return static_cast<void*>(const_cast< JShapeOpDialog*>(this));
    if (!strcmp(_clname, "Ui::ShapeOpDialog"))
        return static_cast< Ui::ShapeOpDialog*>(const_cast< JShapeOpDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JShapeOpDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 17)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 17;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 17)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 17;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
