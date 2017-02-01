/****************************************************************************
** Meta object code from reading C++ file 'MeshFeaturesDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshFeaturesDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshFeaturesDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshFeaturesDialog_t {
    QByteArrayData data[12];
    char stringdata0[196];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshFeaturesDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshFeaturesDialog_t qt_meta_stringdata_JMeshFeaturesDialog = {
    {
        QT_MOC_LITERAL(0, 0, 19), // "JMeshFeaturesDialog"
        QT_MOC_LITERAL(1, 20, 11), // "angleDefect"
        QT_MOC_LITERAL(2, 32, 0), // ""
        QT_MOC_LITERAL(3, 33, 10), // "sharpEdges"
        QT_MOC_LITERAL(4, 44, 20), // "setGaussianCurvature"
        QT_MOC_LITERAL(5, 65, 16), // "setMeanCurvature"
        QT_MOC_LITERAL(6, 82, 15), // "getEigenVectors"
        QT_MOC_LITERAL(7, 98, 18), // "displayEigenVector"
        QT_MOC_LITERAL(8, 117, 22), // "getCurvatureDirections"
        QT_MOC_LITERAL(9, 140, 20), // "openMinKAttribDialog"
        QT_MOC_LITERAL(10, 161, 20), // "openMaxKAttribDialog"
        QT_MOC_LITERAL(11, 182, 13) // "setNodesColor"

    },
    "JMeshFeaturesDialog\0angleDefect\0\0"
    "sharpEdges\0setGaussianCurvature\0"
    "setMeanCurvature\0getEigenVectors\0"
    "displayEigenVector\0getCurvatureDirections\0"
    "openMinKAttribDialog\0openMaxKAttribDialog\0"
    "setNodesColor"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshFeaturesDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    10,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   64,    2, 0x08 /* Private */,
    3,    0,   65,    2, 0x08 /* Private */,
    4,    0,   66,    2, 0x08 /* Private */,
    5,    0,   67,    2, 0x08 /* Private */,
    6,    0,   68,    2, 0x08 /* Private */,
    7,    0,   69,    2, 0x08 /* Private */,
    8,    0,   70,    2, 0x08 /* Private */,
    9,    0,   71,    2, 0x08 /* Private */,
    10,    0,   72,    2, 0x08 /* Private */,
    11,    0,   73,    2, 0x08 /* Private */,

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

    0        // eod
};

void JMeshFeaturesDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshFeaturesDialog *_t = static_cast<JMeshFeaturesDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->angleDefect();
            break;
        case 1:
            _t->sharpEdges();
            break;
        case 2:
            _t->setGaussianCurvature();
            break;
        case 3:
            _t->setMeanCurvature();
            break;
        case 4:
            _t->getEigenVectors();
            break;
        case 5:
            _t->displayEigenVector();
            break;
        case 6:
            _t->getCurvatureDirections();
            break;
        case 7:
            _t->openMinKAttribDialog();
            break;
        case 8:
            _t->openMaxKAttribDialog();
            break;
        case 9:
            _t->setNodesColor();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshFeaturesDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshFeaturesDialog.data,
        qt_meta_data_JMeshFeaturesDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshFeaturesDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshFeaturesDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshFeaturesDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshFeaturesDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshFeaturesDialog"))
        return static_cast< Ui::MeshFeaturesDialog*>(const_cast< JMeshFeaturesDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshFeaturesDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
